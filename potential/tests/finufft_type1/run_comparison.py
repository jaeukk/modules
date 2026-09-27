#!/usr/bin/env python3
"""Build and execute serial, reproducible paired optimization benchmarks."""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
NATIVE = REPO.parent.parent / "control/native-finufft"
DEPS = Path("/home/jaeukk/install")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def build(output, env, entry=None, extra_sources=(), extra_includes=()):
    work = output / "build"
    work.mkdir()
    sources = [entry or HERE / "compare_optimization.cpp"] + [REPO / p for p in (
        "core/GeometryVector.cpp", "core/etc.cpp", "core/RandomGenerator.cpp",
        "core/PeriodicCellList.cpp", "core/SpherePacking.cpp", "core/ConfigurationIO.cpp",
        "potential/Potential.cpp", "potential/StructureOptimization.cpp")] + list(extra_sources)
    common = ["g++", "-std=c++17", "-O2", "-fopenmp", "-DCCPOTENTIAL_USE_FINUFFT",
              "-I" + str(REPO / "core"), "-I" + str(REPO / "potential"),
              "-I" + str(DEPS / "include"), "-I" + str(NATIVE / "install/include"),
              "-I" + str(NATIVE / "source/finufft-2.5.1/include")]
    common.extend("-I" + str(path) for path in extra_includes)
    commands = []

    def compile_one(source):
        obj = work / (source.stem + ".o")
        command = common + ["-c", str(source), "-o", str(obj)]
        with (work / (source.stem + ".log")).open("w") as log:
            result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
        if result.returncode:
            raise RuntimeError("compile failed: " + str(source))
        return obj, command

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        objects = list(pool.map(compile_one, sources))
    commands.extend(command for _, command in objects)
    binary = work / "compare_optimization"
    command = ["g++", "-fopenmp"] + [str(obj) for obj, _ in objects] + [
        "-Wl,--wrap=nlopt_optimize", "-L" + str(DEPS / "lib"),
        "-Wl,-rpath," + str(DEPS / "lib"), "-lnlopt", "-lgsl", "-lgslcblas",
        str(NATIVE / "install/lib/libfinufft.so"),
        "-Wl,-rpath," + str(NATIVE / "install/lib"), "-o", str(binary)]
    commands.append(command)
    save_json(work / "commands.json", commands)
    with (work / "link.log").open("w") as log:
        subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=60)
    save_json(work / "commands.json", commands)
    save_json(work / "source_hashes.json", {str(p): digest(p) for p in sources + [REPO / "potential/CollectiveCoordinatePotential.h"]})
    return binary


def summarize(results):
    rows, speeds, position_rms, micro = [], [], [], []
    all_pass = True
    for data in results:
        initial = data["initial_validation"]
        all_pass &= initial["rho_backend_error"] < 1e-8 and initial["force_backend_error"] < 1e-8
        all_pass &= initial["energy_backend_error"] < 1e-9 * (1 + abs(initial["energy"]))
        for pair in data["pairs"]:
            d, f = pair["direct"], pair["finufft"]
            speed = d["wall_seconds"] / f["wall_seconds"]
            speeds.append(speed)
            position_rms.append(pair["positions"]["aligned_rms"])
            same_input_errors_ok = True
            for run in (d, f):
                v = run["direct_validation"]
                same_input_errors_ok &= v["rho_backend_error"] < 1e-8 and v["force_backend_error"] < 1e-8
                all_pass &= v["ground_state"] and all(0 < p["status"] < 5 for p in run["passes"])
            all_pass &= same_input_errors_ok and f["finufft_evaluations"] == f["rho_evaluations"]
            all_pass &= d["finufft_evaluations"] == 0
            rows.append({"seed": data["seed"], "repeat": pair["repeat"],
                         "direct_seconds": d["wall_seconds"], "finufft_seconds": f["wall_seconds"],
                         "speedup": speed, "direct_energy": d["direct_validation"]["energy"],
                         "finufft_energy_direct_reference": f["direct_validation"]["energy"],
                         "direct_evaluations": sum(p["evaluations"] for p in d["passes"]),
                         "finufft_evaluations": sum(p["evaluations"] for p in f["passes"]),
                         "aligned_rms": pair["positions"]["aligned_rms"]})
        for label, b in [("initial", data["initial_microbench"]),
                         ("direct_endpoint", data["pairs"][0]["direct_endpoint_microbench"]),
                         ("finufft_endpoint", data["pairs"][0]["finufft_endpoint_microbench"])]:
            row = {"seed": data["seed"], "configuration": label}
            for kind in ("rho_seconds", "cached_force_seconds", "energy_force_seconds"):
                direct = statistics.median(b["direct"][kind])
                fast = statistics.median(b["finufft"][kind])
                row[kind] = {"direct_median": direct, "finufft_median": fast, "speedup": direct / fast}
            micro.append(row)
    return {"all_acceptance_checks_pass": bool(all_pass), "optimization_pairs": rows,
            "wall_speedup_median": statistics.median(speeds), "wall_speedup_range": [min(speeds), max(speeds)],
            "aligned_position_rms_range": [min(position_rms), max(position_rms)], "microbenchmarks": micro}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--binary", type=Path, help="Reuse a previously built experiment binary")
    parser.add_argument("--seeds", nargs="+", type=int, default=[20260927, 20260928, 20260929, 20260930, 20261001])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cpus", nargs=2, type=int, help="Optional explicit CPU pair for threading diagnostics")
    parser.add_argument("--omp-proc-bind", choices=["false", "close"], default="false")
    parser.add_argument("--omp-wait-policy", choices=["PASSIVE", "ACTIVE"], default="PASSIVE")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)  # Never overwrite an earlier run.
    available = sorted(os.sched_getaffinity(0))
    # Select distinct physical cores, not neighboring SMT siblings.
    core_ids, physical_cpus = set(), []
    for cpu in available:
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = ((topology / "physical_package_id").read_text().strip(),
                (topology / "core_id").read_text().strip())
        if core not in core_ids:
            core_ids.add(core)
            physical_cpus.append(cpu)
    cpus = args.cpus if args.cpus else physical_cpus[:2]
    if len(cpus) != 2:
        raise RuntimeError("two CPUs required")
    if not set(cpus).issubset(available) or len(set(cpus)) != 2:
        raise RuntimeError("requested CPU pair unavailable")
    os.sched_setaffinity(0, cpus)
    env = dict(os.environ, OMP_NUM_THREADS="2", OMP_THREAD_LIMIT="2", OMP_DYNAMIC="false",
               OMP_PROC_BIND=args.omp_proc_bind, OMP_PLACES="cores", OMP_WAIT_POLICY=args.omp_wait_policy,
               OPENBLAS_NUM_THREADS="2", MKL_NUM_THREADS="2")
    binary = args.binary.resolve() if args.binary else build(output, env)
    settings = {"N": 500, "dimension": 2, "number_density": 1, "box_side": 500**0.5,
                "independent_modes": 299, "chi": 299/998, "seeds": args.seeds, "repeats": args.repeats,
                "optimizer": "unchanged RelaxStructure_NLOPT; NLopt L-BFGS; scales 10000,1,0.01",
                "vector_storage": 100, "max_evaluations_per_pass": 2000,
                "stopping": "repository defaults, no added ftol/xtol; stopval=-MaxEnergy",
                "acceptance": {"direct_energy_max": 1e-11, "direct_max_S": 1e-11, "direct_max_force_component": 1e-7},
                "finufft_tolerance": 1e-12, "cpu_affinity": cpus, "platform": platform.platform(),
                "thread_environment": {k: env[k] for k in env if k.startswith(("OMP_", "OPENBLAS_", "MKL_"))},
                "binary": str(binary), "binary_sha256": digest(binary),
                "header_sha256": digest(REPO / "potential/CollectiveCoordinatePotential.h"),
                "runner_sha256": digest(Path(__file__)),
                "finufft_library_sha256": digest(NATIVE / "install/lib/libfinufft.so"),
                "linked_libraries": subprocess.check_output(["ldd", str(binary)], text=True),
                "compiler": subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
                "position_comparison": "fixed particle labels, periodic minimum-image, best global translation; no permutation matching",
                "mode_selection": "299 shortest half-plane modes; outer shell may be partially selected, as in unit test"}
    save_json(output / "settings.json", settings)
    results = []
    for seed in args.seeds:
        print(f"Running seed {seed}, {args.repeats} paired repetitions", flush=True)
        directory = output / str(seed)
        start = time.monotonic()
        with (output / f"{seed}.json").open("w") as stdout, (output / f"{seed}.stderr.log").open("w") as stderr:
            run = subprocess.run([str(binary), str(directory), str(seed), str(args.repeats)],
                                 env=env, cwd=output, stdout=stdout, stderr=stderr, timeout=600)
        print((output / f"{seed}.stderr.log").read_text(), end="", flush=True)
        if run.returncode not in (0, 2):
            raise RuntimeError(f"seed {seed} failed, exit {run.returncode}; inspect logs")
        result = json.loads((output / f"{seed}.json").read_text())
        result["process_seconds"] = time.monotonic()-start
        result["initial_sha256"] = digest(directory / "initial.tsv")
        results.append(result)
        save_json(output / "summary.json", summarize(results))
    summary = summarize(results)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("optimization_pairs", "microbenchmarks")}, indent=2))
    return 0 if summary["all_acceptance_checks_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
