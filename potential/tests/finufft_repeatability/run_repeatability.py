#!/usr/bin/env python3
"""FINUFFT nthreads=1 repeatability regression.

Does a full fresh build of the current sources, runs API tests with and without the
FINUFFT macro, a 100-repeat fixed-point probe, and 12 serial optimizations (three
seeds, two separate-process repeats per backend). Exits 1 if any required check fails
and 2 on invalid arguments. All artifacts go to a new --output directory.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
TESTS = HERE.parent


def load(name, path):
    """Import a sibling test helper from its file path, without changing sys.path.

    Parameters
    ----------
    name : str
        Module name to register.
    path : Path
        Python file to execute.

    Returns
    -------
    module
        The loaded module.
    """
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


paired = load("finufft_type1_run_comparison", TESTS / "finufft_type1/run_comparison.py")
study = load("five_component_run_study", TESTS / "five_component/run_study.py")
contact = load("contact_reproducibility_run", TESTS / "contact_reproducibility/run.py")

REPO = paired.REPO
APP = REPO / "apps/CCOptimization"
ENTRY = TESTS / "contact_reproducibility/experiment.cpp"
LIBFINUFFT = paired.NATIVE / "lib/libfinufft.so"
SEEDS = [20260928, 20260929, 20260930]
REPEATS = 2
BACKENDS = ["direct", "finufft"]
N, SPECIES, THREADS, MODE_PAIRS = 2000, 5, 4, 80
TARGET_CHI = "0.1"  # The harness selects the complete-shell set nearest this: 80 modes, chi_i = 80/798.
RUN_TIMEOUT = 660
# finufft_observer.cpp wraps the guru calls; ShiftedCCPotential keeps persistent plans.
OBSERVER_WRAP = "-Wl,--wrap=finufft_makeplan,--wrap=finufft_setpts,--wrap=finufft_execute"
RECEIPT = re.compile(r"^FINUFFT_OBSERVER calls=(\d+) null_options=(\d+) nonzero_status=(\d+) nthreads=(\S+)$", re.M)


def save(path, value):
    """Atomically write JSON evidence.

    Parameters
    ----------
    path : Path
        Destination file.
    value : object
        JSON-serializable data.
    """
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


class Checks:
    """Named required checks; the run fails closed if any is false or none was recorded."""

    def __init__(self):
        self.items = []

    def require(self, name, ok, detail=None):
        """Record one required check and print failures immediately.

        Parameters
        ----------
        name : str
            Unique check name.
        ok : bool
            Whether the check passed.
        detail : object, optional
            JSON-serializable context kept with the record.

        Returns
        -------
        bool
            The recorded outcome.
        """
        ok = bool(ok)
        self.items.append(dict(name=name, passed=ok, detail=detail))
        if not ok:
            print(f"FAIL {name}: {detail}", flush=True)
        return ok

    @property
    def passed(self):
        """bool: True only if at least one check was recorded and all passed."""
        return bool(self.items) and all(item["passed"] for item in self.items)

    def failed(self):
        """Return the names of failed checks."""
        return [item["name"] for item in self.items if not item["passed"]]


def observer(text):
    """Parse the single finufft_observer receipt from a stderr log.

    Parameters
    ----------
    text : str
        Captured stderr.

    Returns
    -------
    dict or None
        Call count, null-option and nonzero-status counts, and the nthreads histogram
        (string keys), or None if there is not exactly one receipt.
    """
    found = RECEIPT.findall(text)
    if len(found) != 1:
        return None
    calls, null, bad, histogram = found[0]
    parsed = {} if histogram == "none" else dict(part.split(":") for part in histogram.split(","))
    return dict(calls=int(calls), null_options=int(null), nonzero_status=int(bad),
                nthreads={key: int(value) for key, value in parsed.items()})


def only_one_thread(receipt):
    """Return True if every observed FINUFFT call succeeded and requested nthreads=1."""
    return (receipt is not None and receipt["null_options"] == 0 and receipt["nonzero_status"] == 0
            and receipt["nthreads"] == ({"1": receipt["calls"]} if receipt["calls"] else {}))


def execute(command, stem, env, cwd, timeout):
    """Run one process with stdout/stderr saved next to ``stem``.

    Parameters
    ----------
    command : list of str
        Command line.
    stem : Path
        Log path prefix; ``.stdout`` and ``.stderr`` are appended.
    env : dict
        Process environment.
    cwd : Path
        Working directory.
    timeout : float
        Seconds before the process is killed.

    Returns
    -------
    tuple
        (return code or None on timeout, stdout text, stderr text, wall seconds).
    """
    start = time.monotonic()
    with stem.with_suffix(".stdout").open("w") as out, stem.with_suffix(".stderr").open("w") as err:
        try:
            code = subprocess.run(command, env=env, cwd=cwd, stdout=out, stderr=err, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            code = None
    seconds = time.monotonic() - start
    return code, stem.with_suffix(".stdout").read_text(), stem.with_suffix(".stderr").read_text(), seconds


def tracked_sources():
    """Hash every source, header, and helper that can affect the build or the analysis.

    Returns
    -------
    dict
        Absolute path -> SHA-256.
    """
    files = []
    for directory in (REPO / "core", REPO / "potential"):
        files += sorted(directory.glob("*.h")) + sorted(directory.glob("*.cpp"))
    files += [APP / "MultiStealthyPotential.h", APP / "MultiStealthyPotential.cpp", ENTRY,
              TESTS / "five_component/experiment.cpp", TESTS / "finufft_type1/run_comparison.py",
              TESTS / "five_component/run_study.py", TESTS / "contact_reproducibility/run.py",
              HERE / "finufft_observer.h", HERE / "finufft_observer.cpp", HERE / "api_test.cpp",
              Path(__file__).resolve()]
    return {str(path): paired.digest(path) for path in files}


def build_all(root, env):
    """Fresh full build plus observer relink and both API-test variants.

    The experiment is built by ``run_comparison.build`` from current sources, then relinked
    with the observation-only wrapper. The macro-enabled API test reuses the fresh core
    objects (which do not include CollectiveCoordinatePotential.h) and the fresh
    MultiStealthyPotential object. The macro-disabled API test recompiles only itself and
    MultiStealthyPotential.cpp without the macro and links no FINUFFT library.

    Parameters
    ----------
    root : Path
        Output directory.
    env : dict
        Build environment.

    Returns
    -------
    dict
        Binary paths, build commands, and timings.
    """
    start = time.monotonic()
    plain = paired.build(root, env, ENTRY, [APP / "MultiStealthyPotential.cpp"], [APP])
    commands = json.loads((plain.parent / "commands.json").read_text())
    common = commands[0][:commands[0].index("-c")]
    link = commands[-1]
    if "-DCCPOTENTIAL_USE_FINUFFT" not in common or "-Wl,--wrap=nlopt_optimize" not in link:
        raise RuntimeError("unexpected run_comparison.build command layout")
    objects = {Path(x).stem: x for x in link if x.endswith(".o")}
    extra = root / "extra_build"
    extra.mkdir()
    steps = []

    def step(name, command):
        """Run one extra build command, logging output; raise on failure."""
        command = [str(x) for x in command]
        with (extra / f"{name}.log").open("w") as log:
            code = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=240).returncode
        steps.append(dict(name=name, command=command, returncode=code))
        save(extra / "commands.json", steps)
        if code:
            raise RuntimeError(f"build step {name} failed; see {extra / name}.log")

    observer_obj = extra / "finufft_observer.o"
    step("observer_compile", common + ["-c", HERE / "finufft_observer.cpp", "-o", observer_obj])
    experiment = extra / "experiment_observed"
    relink = list(link)
    relink[relink.index("-o") + 1] = str(experiment)
    at = relink.index("-Wl,--wrap=nlopt_optimize")
    relink[at:at] = [str(observer_obj), OBSERVER_WRAP]
    step("experiment_observed_link", relink)

    core = [objects[name] for name in ("GeometryVector", "etc", "RandomGenerator", "PeriodicCellList",
                                       "SpherePacking", "ConfigurationIO", "Potential")]
    gsl = ["-L" + str(paired.DEPS / "lib"), "-Wl,-rpath," + str(paired.DEPS / "lib"), "-lgsl", "-lgslcblas", "-lm"]
    api_finufft = extra / "api_test_finufft"
    step("api_finufft_compile", common + ["-c", HERE / "api_test.cpp", "-o", extra / "api_test_finufft.o"])
    step("api_finufft_link", ["g++", "-fopenmp", extra / "api_test_finufft.o", *core, objects["MultiStealthyPotential"],
                              observer_obj, OBSERVER_WRAP, *gsl, LIBFINUFFT,
                              "-Wl,-rpath," + str(LIBFINUFFT.parent), "-o", api_finufft])

    # FINUFFT may share the GSL prefix; drop its include path only when it is separate.
    dropped = {"-DCCPOTENTIAL_USE_FINUFFT", "-I" + str(paired.NATIVE / "include")} - {"-I" + str(paired.DEPS / "include")}
    direct_flags = [f for f in common if f not in dropped]
    api_direct = extra / "api_test_direct"
    step("api_direct_compile", direct_flags + ["-c", HERE / "api_test.cpp", "-o", extra / "api_test_direct.o"])
    step("msp_direct_compile", direct_flags + ["-c", APP / "MultiStealthyPotential.cpp",
                                               "-o", extra / "MultiStealthyPotential_direct.o"])
    step("api_direct_link", ["g++", "-fopenmp", extra / "api_test_direct.o", *core,
                             extra / "MultiStealthyPotential_direct.o", *gsl, "-o", api_direct])
    return dict(experiment=experiment, unobserved_experiment=plain, api_finufft=api_finufft, api_direct=api_direct,
                commands=dict(run_comparison_build=commands, extra=steps), seconds=time.monotonic() - start)


def api_tests(root, binaries, env, checks):
    """Run both API-test variants and the no-FINUFFT dependency check.

    Returns
    -------
    dict
        Parsed receipts per variant.
    """
    directory = root / "api"
    directory.mkdir()
    results = {}
    for name, macro in (("api_finufft", True), ("api_direct", False)):
        code, out, err, seconds = execute([str(binaries[name])], directory / name, env, root, 120)
        try:
            data = json.loads(out.strip().splitlines()[-1])
        except (IndexError, ValueError):
            data = {}
        receipt = observer(err)
        results[name] = dict(returncode=code, result=data, observer=receipt, seconds=seconds)
        checks.require(f"{name}_passed", code == 0 and data.get("passed") is True and data.get("finufft_macro") is macro,
                       dict(returncode=code, result=data))
        if macro:
            checks.require(f"{name}_observer", receipt is not None and receipt["calls"] > 0
                           and receipt["calls"] == data.get("observed_finufft_calls")
                           and receipt["null_options"] == receipt["nonzero_status"] == 0, receipt)
        else:
            checks.require(f"{name}_no_observer_receipt", receipt is None, receipt)
    linked = subprocess.run(["ldd", str(binaries["api_direct"])], capture_output=True, text=True).stdout
    results["api_direct_ldd"] = linked
    checks.require("api_direct_has_no_libfinufft", "finufft" not in linked, linked)
    save(directory / "results.json", results)
    return results


def probe(root, experiment, env, checks):
    """Run the contact-reproducibility probe and require identical values over 100 repeats.

    Returns
    -------
    dict
        Parsed probe output and observer receipt.
    """
    code, out, err, seconds = execute([str(experiment), "probe"], root / "probe", env, root, 180)
    try:
        data = json.loads(out)
    except ValueError:
        data = {}
    receipt = observer(err)
    backends = {b.get("backend"): b for b in data.get("backends", [])}
    checks.require("probe_passed", code == 0 and data.get("passed") is True
                   and (data.get("N"), data.get("species"), data.get("threads")) == (N, SPECIES, THREADS),
                   dict(returncode=code, stderr_tail=err[-2000:]))
    for backend in BACKENDS:
        values = backends.get(backend, {})
        checks.require(f"probe_{backend}_100_repeats_identical",
                       values.get("energy_range") == 0.0 and values.get("force_range") == 0.0, values)
    checks.require("probe_finufft_observed_nthreads_1", only_one_thread(receipt) and receipt["calls"] > 0, receipt)
    result = dict(returncode=code, seconds=seconds, result=data, observer=receipt)
    save(root / "probe.json", result)
    return result


def optimize(root, experiment, env, seed, backend, repeat, checks):
    """Run one optimization in its own process and apply per-run required checks.

    Returns
    -------
    dict
        Receipt with command, timings, harness output, observer receipt, and file hashes.
    """
    name = f"{seed}_{backend}_{repeat}"
    dest = root / "runs" / name
    command = [str(experiment), str(dest), str(seed), TARGET_CHI, backend]
    print(f"run {name}", flush=True)
    code, out, err, seconds = execute(command, root / "runs" / name, env, root, RUN_TIMEOUT)
    row = dict(name=name, seed=seed, selected_backend=backend, repeat=repeat, command=command,
               returncode=code, process_seconds=seconds, observer=observer(err))
    try:
        row["result"] = json.loads(out)
    except ValueError:
        row["result"] = None
    for file in ("initial.tsv", "modes.tsv", "final.tsv", "direct_structure.tsv"):
        if (dest / file).exists():
            row[file.replace(".", "_") + "_sha256"] = paired.digest(dest / file)
    result = row["result"] or {}
    reference = result.get("direct_reference", {})
    passes = result.get("passes", [])
    rho, fast = result.get("rho_calls"), result.get("finufft_calls")
    row["actual_backend"] = ("finufft" if rho and fast == rho else "direct" if fast == 0
                             else "mixed" if fast else "unknown")
    row["optimizer_statuses"] = [p["status"] for p in passes]
    row["evaluations"] = sum(p["evaluations"] for p in passes)
    row["optimizer_codes_success"] = bool(passes) and all(0 < p["status"] < 5 for p in passes)
    checks.require(f"{name}_exit_zero", code == 0 and row["result"] is not None,
                   dict(returncode=code, stderr_tail=err[-2000:]))
    checks.require(f"{name}_settings", result.get("backend") == backend and result.get("N") == N
                   and result.get("threads") == THREADS and result.get("mode_pairs") == MODE_PAIRS
                   and abs(result.get("actual_chi", 0) - MODE_PAIRS / (2 * (N // SPECIES - 1))) < 1e-15,
                   {k: result.get(k) for k in ("backend", "N", "threads", "mode_pairs", "actual_chi")})
    checks.require(f"{name}_direct_reference_accepted", reference.get("ground_state") is True, reference)
    checks.require(f"{name}_no_time_budget_truncation", bool(passes) and all(p["status"] != 6 for p in passes)
                   and result.get("optimization_seconds", 600) < 600, row["optimizer_statuses"])
    if backend == "finufft":
        checks.require(f"{name}_finufft_used_without_fallback", row["actual_backend"] == "finufft",
                       dict(rho_calls=rho, finufft_calls=fast))
        checks.require(f"{name}_observed_nthreads_1", only_one_thread(row["observer"])
                       and row["observer"]["calls"] >= (fast or 1), row["observer"])
    else:
        checks.require(f"{name}_direct_made_no_finufft_calls", fast == 0 and row["observer"] is not None
                       and row["observer"]["calls"] == 0, dict(finufft_calls=fast, observer=row["observer"]))
    save(root / "runs" / f"{name}.json", row)
    return row


def compare_seed(root, seed, rows, checks):
    """Require identical inputs and same-backend repeats; describe cross-backend displacements.

    Returns
    -------
    dict
        Identity flags and displacement statistics for one seed.
    """
    runs = root / "runs"
    group = {(r["selected_backend"], r["repeat"]): r for r in rows if r["seed"] == seed}
    record = dict(seed=seed)
    for file in ("initial_tsv_sha256", "modes_tsv_sha256"):
        values = {group[k].get(file) for k in group}
        checks.require(f"{seed}_identical_{file[:-7]}", len(group) == 4 and len(values) == 1 and None not in values, sorted(map(str, values)))
    record["same_backend"] = {}
    for backend in BACKENDS:
        a, b = group.get((backend, 0)), group.get((backend, 1))
        ra, rb = (a or {}).get("result") or {}, (b or {}).get("result") or {}
        final_a, final_b = runs / f"{seed}_{backend}_0/final.tsv", runs / f"{seed}_{backend}_1/final.tsv"
        identical = final_a.exists() and final_b.exists() and final_a.read_bytes() == final_b.read_bytes()
        checks.require(f"{seed}_{backend}_final_byte_identical", identical)
        fields = ("passes", "rho_calls", "finufft_calls")
        checks.require(f"{seed}_{backend}_pass_records_identical",
                       bool(ra) and all(ra.get(f) == rb.get(f) for f in fields), {f: [ra.get(f), rb.get(f)] for f in fields})
        checks.require(f"{seed}_{backend}_observer_receipts_identical",
                       a is not None and b is not None and a["observer"] == b["observer"],
                       [(a or {}).get("observer"), (b or {}).get("observer")])
        entry = dict(final_byte_identical=identical)
        if final_a.exists() and final_b.exists():
            try:
                entry["displacement"] = contact.compare(final_a, final_b, N)
            except AssertionError as error:
                entry["displacement_error"] = f"label/species/radius mismatch {error}"
        record["same_backend"][backend] = entry
    record["cross_backend_descriptive"] = []
    for repeat in range(REPEATS):
        final_d, final_f = runs / f"{seed}_direct_{repeat}/final.tsv", runs / f"{seed}_finufft_{repeat}/final.tsv"
        if final_d.exists() and final_f.exists():
            try:
                record["cross_backend_descriptive"].append(dict(repeat=repeat, **contact.compare(final_d, final_f, N)))
            except AssertionError as error:
                record["cross_backend_descriptive"].append(dict(repeat=repeat, error=f"label/species/radius mismatch {error}"))
    return record


def report(root, summary):
    """Write a short human-readable REPORT.md from summary.json content."""
    lines = ["# FINUFFT nthreads=1 repeatability regression", "",
             f"Result: **{'PASS' if summary['passed'] else 'FAIL'}** ({len(summary['checks'])} required checks, "
             f"{len(summary['failed_checks'])} failed).", "",
             "N=2000, five equal components, phi_total=0.2, chi_i=80/798, OMP_NUM_THREADS=OMP_THREAD_LIMIT=4, "
             f"passive waiting, CPU affinity {summary.get('cpu_affinity')}. Production defaults only: the FINUFFT "
             "option is never set by the harness, so these runs exercise the default nthreads=1.", "",
             "Finite regression evidence on this machine, build, and settings only; not a general guarantee. "
             "Cross-backend displacements are descriptive, not acceptance criteria.", ""]
    if summary["failed_checks"]:
        lines += ["Failed checks: " + ", ".join(summary["failed_checks"]), ""]
    lines += ["| Run | Actual backend | Optimization s | Process s | Evaluations | Status codes | Accepted | FINUFFT calls (harness/observed) |",
              "|---|---|---:|---:|---:|---|---|---|"]
    for row in summary.get("runs", []):
        result = row.get("result") or {}
        observed = (row.get("observer") or {}).get("calls")
        lines.append(f"| {row['name']} | {row['actual_backend']} | {result.get('optimization_seconds', float('nan')):.3f} | "
                     f"{row['process_seconds']:.3f} | {row['evaluations']} | {','.join(map(str, row['optimizer_statuses']))} | "
                     f"{result.get('direct_reference', {}).get('ground_state')} | {result.get('finufft_calls')}/{observed} |")
    lines += ["", "| Seed | Direct repeats identical | FINUFFT repeats identical | Cross-backend mean displacement r0/r1 |",
              "|---|---|---|---|"]
    for record in summary.get("comparisons", []):
        cross = "/".join(f"{c['mean']:.3e}" if "mean" in c else "n/a" for c in record["cross_backend_descriptive"])
        same = record["same_backend"]
        lines.append(f"| {record['seed']} | {same['direct']['final_byte_identical']} | {same['finufft']['final_byte_identical']} | {cross} |")
    for backend, times in summary.get("timing", {}).items():
        lines.append("")
        lines.append(f"{backend}: median optimization {times['median_optimization_seconds']:.3f} s, "
                     f"range {times['optimization_seconds_range'][0]:.3f}-{times['optimization_seconds_range'][1]:.3f} s.")
    (root / "REPORT.md").write_text("\n".join(lines) + "\n")


def main():
    """Parse arguments, run the whole regression serially, and return the exit code."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, required=True, help="new absolute directory outside the source tree")
    args = parser.parse_args()
    if not args.output.is_absolute():
        parser.error("--output must be an absolute path")
    root = args.output.resolve()
    source_root = REPO.parent.resolve()
    if root == source_root or source_root in root.parents:
        parser.error(f"--output must be outside the source tree {source_root}")
    if root.exists():
        parser.error(f"refusing to overwrite existing {root}")
    root.mkdir(parents=True, exist_ok=False)
    checks = Checks()
    summary = dict(started=time.time(), passed=False, checks=checks.items, failed_checks=[])
    try:
        env, cpus = study.environment()
        summary.update(cpu_affinity=cpus, thread_environment={k: v for k, v in env.items() if k.startswith(("OMP_", "OPENBLAS_", "MKL_"))})
        checks.require("thread_environment", env["OMP_NUM_THREADS"] == env["OMP_THREAD_LIMIT"] == "4"
                       and env["OMP_WAIT_POLICY"] == "PASSIVE" and len(set(cpus)) == 4, summary["thread_environment"])
        sources_before = tracked_sources()
        built = build_all(root, env)
        binaries = {k: built[k] for k in ("experiment", "unobserved_experiment", "api_finufft", "api_direct")}
        manifest = dict(platform=platform.platform(), python=sys.version,
                        compiler=subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
                        cpu_affinity=cpus, thread_environment=summary["thread_environment"],
                        build_seconds=built["seconds"], commands=built["commands"], sources=sources_before,
                        binaries={k: dict(path=str(v), sha256=paired.digest(v)) for k, v in binaries.items()},
                        libfinufft=dict(path=str(LIBFINUFFT), sha256=paired.digest(LIBFINUFFT)),
                        ldd={k: subprocess.run(["ldd", str(v)], capture_output=True, text=True).stdout for k, v in binaries.items()},
                        settings=dict(N=N, species=SPECIES, threads=THREADS, packing_fraction=0.2, target_chi=TARGET_CHI,
                                      mode_pairs=MODE_PAIRS, seeds=SEEDS, repeats=REPEATS, run_timeout_seconds=RUN_TIMEOUT))
        save(root / "manifest.json", manifest)
        summary["api"] = api_tests(root, binaries, env, checks)
        summary["probe"] = probe(root, built["experiment"], env, checks)
        (root / "runs").mkdir()
        rows = []
        for seed in SEEDS:
            for repeat in range(REPEATS):
                for backend in (BACKENDS if (seed + repeat) % 2 == 0 else BACKENDS[::-1]):
                    rows.append(optimize(root, built["experiment"], env, seed, backend, repeat, checks))
                    summary["runs"] = rows
                    save(root / "summary.json", summary)
        summary["comparisons"] = [compare_seed(root, seed, rows, checks) for seed in SEEDS]
        summary["timing"] = {}
        for backend in BACKENDS:
            times = [r["result"]["optimization_seconds"] for r in rows if r["selected_backend"] == backend and r["result"]]
            if times:
                summary["timing"][backend] = dict(median_optimization_seconds=statistics.median(times),
                                                  optimization_seconds_range=[min(times), max(times)],
                                                  process_seconds=[r["process_seconds"] for r in rows if r["selected_backend"] == backend])
        sources_after = tracked_sources()
        checks.require("sources_unchanged_during_run", sources_after == sources_before,
                       sorted(k for k in set(sources_before) | set(sources_after) if sources_before.get(k) != sources_after.get(k)))
    except Exception as error:  # Any build or harness error fails closed with a recorded reason.
        checks.require("driver_completed", False, "".join(traceback.format_exception(type(error), error, error.__traceback__)))
    summary.update(finished=time.time(), passed=checks.passed, failed_checks=checks.failed())
    save(root / "summary.json", summary)
    report(root, summary)
    print(json.dumps(dict(passed=summary["passed"], failed_checks=summary["failed_checks"], output=str(root)), indent=2))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
