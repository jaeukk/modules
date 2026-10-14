"""Benchmark matrix for bench_finufft.out: per-evaluation timings, optimization wall time and
same-backend repeatability (two separate processes per backend must give identical coordinates).

Usage: python3 tests/bench_finufft.py [--threads 4] [--seed 1] [--quick] [--json OUT]
Build first: make USE_FINUFFT=1 bench_finufft.out. Prints a markdown table; --json writes raw rows
(put it outside the repository).
"""
import argparse
import json
import os
import pathlib
import subprocess

EXE = pathlib.Path(__file__).resolve().parents[1] / "bench_finufft.out"
CASES = [(2, 0.1), (2, 0.4), (3, 0.1), (3, 0.4)]


def run(*args, threads):
    """Run the benchmark executable and return its JSON rows."""
    env = dict(os.environ, OMP_PROC_BIND="close", OMP_PLACES="cores")
    out = subprocess.run([str(EXE), *map(str, args)], env=env, check=True,
                         capture_output=True, text=True).stdout
    return [json.loads(line) for line in out.splitlines() if line.startswith("{")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--quick", action="store_true", help="2D cases only, no optimization runs")
    ap.add_argument("--json", type=pathlib.Path)
    a = ap.parse_args()
    rows = []
    print("| d | chi_s | modes | backend | eval ms | E rel diff | F rel diff | opt s (r0/r1) | calls | ms/call | final E | repeatable |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for dim, chi in CASES[:2] if a.quick else CASES:
        evals = run("eval", dim, 2000, 5, chi, a.threads, a.seed, threads=a.threads)
        rows += evals
        for ev in evals:
            opt = "" if a.quick else [run("opt", dim, 2000, 5, chi, a.threads, a.seed, ev["backend"],
                                          threads=a.threads)[0] for _ in range(2)]
            rows += opt or []
            o = opt[0] if opt else {}
            cells = [dim, chi, ev["modes"], ev["backend"], f"{ev['median_total_ms']:.3f}",
                     f"{ev['energy_rel_diff']:.1e}", f"{ev['force_rel_diff']:.1e}",
                     f"{opt[0]['seconds']:.2f} / {opt[1]['seconds']:.2f}" if opt else "",
                     o.get("energy_calls", ""), f"{o['ms_per_call']:.3f}" if opt else "",
                     f"{o['final_energy_direct']:.2e}" if opt else "",
                     ("yes" if opt[0]["hash"] == opt[1]["hash"] else "NO") if opt else ""]
            print("| " + " | ".join(map(str, cells)) + " |", flush=True)
    if a.json:
        a.json.write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
