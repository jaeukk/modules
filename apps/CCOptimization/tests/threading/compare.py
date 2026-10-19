#!/usr/bin/env python3
"""Compare Options A and B on the same four CPUs; keep outputs outside datasets."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess


def choose_cpus(count):
    """Select available CPUs, preferring distinct exposed cores.

    Parameters
    ----------
    count : int
        Number of logical CPUs required by the comparison.

    Returns
    -------
    list of int
        CPU IDs available in the current process affinity mask.

    Raises
    ------
    RuntimeError
        If fewer than ``count`` CPUs are available.
    """
    available = sorted(os.sched_getaffinity(0))
    selected, cores = [], set()
    for cpu in available:
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        core = tuple((topology / name).read_text().strip() for name in ("physical_package_id", "core_id"))
        if core not in cores:
            selected.append(cpu)
            cores.add(core)
        if len(selected) == count:
            return selected
    selected.extend(cpu for cpu in available if cpu not in selected)
    if len(selected) < count:
        raise RuntimeError(f"Need {count} available CPUs, found {len(selected)}")
    return selected[:count]


def main():
    """Run matched evaluator and bounded NLopt comparisons.

    Returns
    -------
    None
        Results and provenance are written to the requested new directory.

    Raises
    ------
    RuntimeError
        If a child fails or produces unexpected diagnostic output.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="New directory for metadata, logs, and measurements")
    parser.add_argument("--eval-calls", type=int, default=400)
    parser.add_argument("--opt-maxsteps", type=int, default=1000, help="Evaluation cap per NLopt scaling phase")
    parser.add_argument("--history", type=int, default=40, help="L-BFGS correction pairs; use 100 to reproduce the original A/B comparison")
    args = parser.parse_args()
    if args.eval_calls < 1 or args.opt_maxsteps < 1 or args.history < 1:
        parser.error("Evaluation counts and history must be positive")
    binary = Path(__file__).resolve().parents[2] / "bench_threading.out"
    if not binary.is_file():
        parser.error("Build bench_threading.out with USE_FINUFFT=1 first")
    cpus = choose_cpus(4)
    args.output.mkdir(parents=True, exist_ok=False)
    os.sched_setaffinity(0, cpus)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="4", OMP_THREAD_LIMIT="4", OMP_DYNAMIC="false",
               OMP_WAIT_POLICY="PASSIVE", OMP_PROC_BIND="false", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    repo = Path(__file__).resolve().parents[4]
    metadata = dict(commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
                    dirty=bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=repo, text=True)),
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    cpu_model=next(line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines() if line.startswith("model name")),
                    platform=platform.platform(), affinity=cpus, species=3, dimension=2,
                    particles_per_species=[300, 600, 1000], seeds=[17, 18, 19], threads=4,
                    phi=.70, chi_per_species=.30, contact_strength=100, history=args.history,
                    finufft_tolerance=1e-12, eval_calls=args.eval_calls,
                    opt_maxsteps_per_phase=args.opt_maxsteps,
                    environment={k: env[k] for k in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OMP_DYNAMIC", "OMP_WAIT_POLICY", "OMP_PROC_BIND")},
                    option_A="Concurrent species, internal FINUFFT threads=1",
                    option_B="Sequential species, internal FINUFFT threads=4, extra contact worker time-shares affinity")
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    results = []
    # Interleave policies and rotate size order to reduce order/thermal bias.
    for mode in ("eval", "opt"):
        repetitions = range(3) if mode == "eval" else range(1)
        for seed in metadata["seeds"]:
            for repetition in repetitions:
                sizes = metadata["particles_per_species"]
                offset = (seed + repetition) % len(sizes)
                for size in sizes[offset:] + sizes[:offset]:
                    options = ("A", "B") if (seed + repetition) % 2 else ("B", "A")
                    for option in options:
                        name = f"{mode}-n{size}-{option}-seed{seed}-r{repetition}"
                        folder = args.output / name
                        folder.mkdir()
                        iterations = args.eval_calls if mode == "eval" else args.opt_maxsteps
                        command = [str(binary), mode, option, str(size), "4", str(seed), str(iterations), str(repetition), str(args.history)]
                        (folder / "command.json").write_text(json.dumps(command) + "\n")
                        process = subprocess.run(command, cwd=folder, env=env, capture_output=True, text=True, timeout=300)
                        (folder / "stdout.log").write_text(process.stdout)
                        (folder / "stderr.log").write_text(process.stderr)
                        if process.returncode != 0:
                            raise RuntimeError(f"{name} failed: {process.stderr.strip()}")
                        if process.stderr.strip():
                            raise RuntimeError(f"{name} produced stderr; inspect {folder}")
                        result = json.loads(process.stdout)
                        results.append(result)
                        (args.output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
                        print(json.dumps(result), flush=True)
    summary = []
    for mode in ("eval", "opt"):
        for size in metadata["particles_per_species"]:
            for option in ("A", "B"):
                group = [x for x in results if (x["mode"], x["particles_per_species"], x["option"]) == (mode, size, option)]
                summary.append(dict(mode=mode, particles_per_species=size, option=option,
                                    runs=len(group), **{key: statistics.median(x[key] for x in group)
                                    for key in ("wall_seconds", "cpu_seconds", "active_cpus", "utilization_percent", "evaluations")}))
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Completed {len(results)} runs; results: {args.output}", flush=True)


if __name__ == "__main__":
    main()
