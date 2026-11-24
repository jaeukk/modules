#!/usr/bin/env python3
"""Compare time per optimizer step of GetMultiCCO variants on the same initial states.

Runs ``bench_step.out`` once per (mode, variant, size, seed, repetition) on four pinned CPUs,
interleaving the variants, and writes raw results, provenance and a summary to a new
directory outside the repository. Results of earlier rounds (the controls) can be merged
into the summary with ``--reuse`` instead of being measured again.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
from pathlib import Path

from compare import choose_cpus

THREADS = 4
ENVIRONMENT = dict(OMP_NUM_THREADS=str(THREADS), OMP_THREAD_LIMIT=str(THREADS), OMP_DYNAMIC="false",
                   OMP_WAIT_POLICY="PASSIVE", OMP_PROC_BIND="false", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")


def counts(variant: str, size: int, mode: str, args: argparse.Namespace) -> int:
    """Return the number of calls (eval) or evaluations per scaling phase (opt) of one run.

    Direct sums cost O(N M) per evaluation, so their runs are shortened; time per step
    does not depend on the length of a run.

    Parameters
    ----------
    variant : str
        Variant name, e.g. ``"A+grid"``.
    size : int
        Particles per species.
    mode : str
        ``"eval"`` or ``"opt"``.
    args : argparse.Namespace
        Parsed command line (``eval_calls``, ``opt_maxsteps``).

    Returns
    -------
    int
        Count passed to ``bench_step.out``.
    """
    full = args.eval_calls if mode == "eval" else args.opt_maxsteps
    if "direct" not in variant.split("+"):
        return full
    scale = max(1.0, (size / 300.0) ** 2)
    return max(12 if mode == "eval" else 20, int(full / scale))


def summarize(results: list[dict], baseline: str) -> list[dict]:
    """Group runs by mode, size and variant and reduce each group to medians.

    Parameters
    ----------
    results : list of dict
        JSON records of ``bench_step.out``, each with an added ``round`` label.
    baseline : str
        Variant whose median time per step defines ``speedup`` within a mode and size.

    Returns
    -------
    list of dict
        One row per group: run count, median, minimum and maximum time per step [ms],
        median time inside and outside the potential [ms], CPU utilization [%],
        evaluations, median wall time of a run [s], number of runs whose final energy by
        direct sums is below 1e-14 (accepted ground states), speedup against the baseline,
        and the rounds the runs came from.
    """
    rows = []
    keys = sorted({(r["mode"], r["particles_per_species"], r["variant"]) for r in results})
    for mode, size, variant in keys:
        group = [r for r in results if (r["mode"], r["particles_per_species"], r["variant"]) == (mode, size, variant)]
        rows.append(dict(mode=mode, particles_per_species=size, variant=variant, runs=len(group),
                         ms_per_step=statistics.median(r["ms_per_step"] for r in group),
                         ms_min=min(r["ms_per_step"] for r in group), ms_max=max(r["ms_per_step"] for r in group),
                         potential_ms=statistics.median(r["potential_ms_per_step"] for r in group),
                         optimizer_ms=statistics.median(r["optimizer_ms_per_step"] for r in group),
                         utilization_percent=statistics.median(r["utilization_percent"] for r in group),
                         evaluations=statistics.median(r["evaluations"] for r in group),
                         wall_seconds=statistics.median(r["wall_seconds"] for r in group),
                         accepted=sum(bool(r.get("accepted")) for r in group),
                         rounds=sorted({r["round"] for r in group})))
    for row in rows:
        base = [b for b in rows if (b["mode"], b["particles_per_species"], b["variant"]) ==
                (row["mode"], row["particles_per_species"], baseline)]
        row["speedup"] = base[0]["ms_per_step"] / row["ms_per_step"] if base else None
    return rows


def table(rows: list[dict]) -> str:
    """Format summary rows as a Markdown table.

    Parameters
    ----------
    rows : list of dict
        Output of :func:`summarize`.

    Returns
    -------
    str
        Markdown text, one table per mode.
    """
    lines = []
    for mode, title in (("opt", "Bounded optimization: wall time per objective evaluation"),
                        ("eval", "Warm evaluator: wall time per SetConfiguration + EnergyAndAllForce")):
        lines += [f"### {title}", "",
                  "| N/species | Variant | Runs | ms/step (median) | min-max | In potential | Outside | CPU % | Evaluations | Wall s | Ground states | Speedup | Rounds |",
                  "|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---|"]
        for row in (r for r in rows if r["mode"] == mode):
            speedup = "" if row["speedup"] is None else f"{row['speedup']:.2f}"
            lines.append(f"| {row['particles_per_species']} | {row['variant']} | {row['runs']} | {row['ms_per_step']:.3f} | "
                         f"{row['ms_min']:.3f}-{row['ms_max']:.3f} | {row['potential_ms']:.3f} | {row['optimizer_ms']:.3f} | "
                         f"{row['utilization_percent']:.0f} | {row['evaluations']:.0f} | {row['wall_seconds']:.2f} | "
                         f"{row['accepted']}/{row['runs']} | {speedup} | {','.join(row['rounds'])} |")
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    """Run the requested variants and write results, provenance and the summary.

    Raises
    ------
    RuntimeError
        If a child process fails or writes to standard error.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", required=True, type=Path, help="New directory for metadata, logs and measurements")
    parser.add_argument("--round", required=True, help="Label of this round, recorded with every run")
    parser.add_argument("--variants", required=True, help="Comma-separated variants, e.g. direct,A,B,A+grid,A+fast")
    parser.add_argument("--sizes", default="300,1000,4000", help="Particles per species, comma-separated")
    parser.add_argument("--seeds", default="17,18,19")
    parser.add_argument("--eval-calls", type=int, default=300)
    parser.add_argument("--opt-maxsteps", type=int, default=300, help="Evaluation cap per NLopt scaling phase")
    parser.add_argument("--eval-reps", type=int, default=3)
    parser.add_argument("--opt-reps", type=int, default=2)
    parser.add_argument("--history", type=int, default=40)
    parser.add_argument("--contact", type=float, default=100.0, help="Relative contact strength of the fixture")
    parser.add_argument("--phi", type=float, default=0.70, help="Packing fraction of the fixture")
    parser.add_argument("--modes", default="eval,opt", help="Comma-separated subset of eval,opt")
    parser.add_argument("--baseline", default="A", help="Variant against which speedups are reported")
    parser.add_argument("--reuse", type=Path, nargs="*", default=[], help="results.json files of earlier rounds to merge into the summary")
    args = parser.parse_args()
    variants = args.variants.split(",")
    sizes = [int(x) for x in args.sizes.split(",")]
    seeds = [int(x) for x in args.seeds.split(",")]
    binary = Path(__file__).resolve().parents[2] / "bench_step.out"
    if not binary.is_file():
        parser.error("Build bench_step.out with USE_FINUFFT=1 first")
    cpus = choose_cpus(THREADS)
    args.output.mkdir(parents=True, exist_ok=False)
    os.sched_setaffinity(0, cpus)
    env = os.environ.copy()
    env.update(ENVIRONMENT)
    repo = Path(__file__).resolve().parents[4]
    metadata = dict(round=args.round, commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
                    dirty=bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=repo, text=True)),
                    diff_sha256=hashlib.sha256(subprocess.check_output(["git", "diff", "HEAD"], cwd=repo)).hexdigest(),
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    cpu_model=next(line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines() if line.startswith("model name")),
                    platform=platform.platform(), affinity=cpus, load_average=os.getloadavg(), species=3, dimension=2,
                    particles_per_species=sizes, seeds=seeds, threads=THREADS, phi=args.phi, chi_per_species=.30,
                    contact_strength=args.contact, history=args.history, finufft_tolerance=1e-12, variants=variants,
                    eval_calls=args.eval_calls, opt_maxsteps_per_phase=args.opt_maxsteps, eval_reps=args.eval_reps,
                    opt_reps=args.opt_reps, environment=ENVIRONMENT, reused=[str(p) for p in args.reuse])
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    results = []
    # Interleave the variants and rotate their order and the size order to spread drift over all of them.
    for mode in args.modes.split(","):
        for seed in seeds:
            for repetition in range(args.eval_reps if mode == "eval" else args.opt_reps):
                shift = (seed + repetition) % len(sizes)
                for size in sizes[shift:] + sizes[:shift]:
                    turn = (seed + repetition + size) % len(variants)
                    for variant in variants[turn:] + variants[:turn]:
                        name = f"{mode}-n{size}-{variant}-seed{seed}-r{repetition}"
                        folder = args.output / name
                        folder.mkdir()
                        command = [str(binary), mode, variant, str(size), str(THREADS), str(seed),
                                   str(counts(variant, size, mode, args)), str(repetition), str(args.history),
                                   str(args.contact), str(args.phi)]
                        # "+active" is a setting, not code: the OpenMP threads spin between parallel regions instead of sleeping.
                        child = dict(env, OMP_WAIT_POLICY="ACTIVE") if "active" in variant.split("+") else env
                        process = subprocess.run(command, cwd=folder, env=child, capture_output=True, text=True, timeout=1800)
                        (folder / "stdout.log").write_text(process.stdout)
                        if process.returncode != 0 or process.stderr.strip():
                            (folder / "stderr.log").write_text(process.stderr)
                            raise RuntimeError(f"{name} failed: {process.stderr.strip()}")
                        result = json.loads(process.stdout)
                        result["round"] = args.round
                        results.append(result)
                        (args.output / "results.json").write_text(json.dumps(results, indent=1) + "\n")
                        print(f"{name}: {result['ms_per_step']:.3f} ms/step", flush=True)
    merged = list(results)
    for path in args.reuse:
        merged += json.loads(path.read_text())
    rows = summarize(merged, args.baseline)
    (args.output / "summary.json").write_text(json.dumps(rows, indent=1) + "\n")
    text = table(rows)
    (args.output / "summary.md").write_text(text + "\n")
    print(text)
    print(f"Completed {len(results)} runs; results: {args.output}", flush=True)


if __name__ == "__main__":
    main()
