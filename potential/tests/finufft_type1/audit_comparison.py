#!/usr/bin/env python3
"""Independent scalar direct-sum audit of saved experiment data; no C++/FINUFFT calls."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics


def rows(path, cast=float):
    with path.open() as file:
        return [tuple(map(cast, line.split())) for line in file if line.strip()]


def direct_reference(positions, modes, side):
    n = len(positions)
    volume = side * side
    forces = [[0.0, 0.0] for _ in positions]
    structure = []
    for mx, my in modes:
        phase = [2 * math.pi * (mx*x + my*y) for x, y in positions]
        cosine, sine = [math.cos(x) for x in phase], [math.sin(x) for x in phase]
        real, imag = math.fsum(cosine), math.fsum(sine)
        structure.append((real*real + imag*imag)/n)
        for j, (co, si) in enumerate(zip(cosine, sine)):
            coefficient = 2/volume * (real*si - imag*co) * (2*math.pi/side)
            forces[j][0] += coefficient*mx
            forces[j][1] += coefficient*my
    return {"energy": math.fsum(structure)*n/volume, "max_S": max(structure),
            "max_force_component": max(abs(x) for force in forces for x in force)}


def periodic(x):
    return x-math.floor(x+0.5)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    root = args.directory.resolve()
    settings = json.loads((root/"settings.json").read_text())
    summary = json.loads((root/"summary.json").read_text())
    side = settings["box_side"]
    audits, table, all_pairs, all_micro = [], [], [], []
    for seed in settings["seeds"]:
        data = json.loads((root/f"{seed}.json").read_text())
        directory = root/str(seed)
        modes = rows(directory/"modes.tsv", int)
        assert len(modes) == 299 and len(set(modes)) == 299
        assert not any((-x, -y) in modes for x, y in modes)
        for pair in data["pairs"]:
            rep = pair["repeat"]
            endpoints = {}
            for backend in ("direct", "finufft"):
                positions = rows(directory/f"{backend}_{rep}.tsv")
                assert len(positions) == 500 and all(len(p) == 2 for p in positions)
                assert all(math.isfinite(v) and 0 <= v < 1 for p in positions for v in p)
                endpoints[backend] = positions
                reference = direct_reference(positions, modes, side)
                recorded = pair[backend]["direct_validation"]
                assert reference["energy"] <= 1e-11 and reference["max_S"] <= 1e-11
                assert reference["max_force_component"] <= 1e-7
                assert abs(reference["energy"]-recorded["energy"]) <= 1e-18
                assert abs(reference["max_force_component"]-recorded["max_force_component"]) <= 1e-9
                audits.append({"seed": seed, "repeat": rep, "backend": backend, **reference,
                               "energy_difference_from_cpp": reference["energy"]-recorded["energy"],
                               "max_force_difference_from_cpp": reference["max_force_component"]-recorded["max_force_component"]})
            shift = pair["positions"]["translation_cartesian"]
            distances = [side*math.hypot(*(periodic(b[d]-a[d]-shift[d]/side) for d in range(2)))
                         for a, b in zip(endpoints["direct"], endpoints["finufft"])]
            rms = math.sqrt(math.fsum(x*x for x in distances)/500)
            assert abs(rms-pair["positions"]["aligned_rms"]) < 1e-14
            assert abs(max(distances)-pair["positions"]["aligned_max"]) < 1e-14
            assert rms <= pair["positions"]["raw_periodic_rms"] + 1e-14
            all_pairs.append(pair)
        pairs = data["pairs"]
        table.append({"seed": seed,
                      "direct_seconds": statistics.median(p["direct"]["wall_seconds"] for p in pairs),
                      "finufft_seconds": statistics.median(p["finufft"]["wall_seconds"] for p in pairs),
                      "direct_energy_max": max(p["direct"]["direct_validation"]["energy"] for p in pairs),
                      "finufft_energy_max": max(p["finufft"]["direct_validation"]["energy"] for p in pairs),
                      "aligned_rms_max": max(p["positions"]["aligned_rms"] for p in pairs),
                      "direct_evaluations": sorted(set(sum(s["evaluations"] for s in p["direct"]["passes"]) for p in pairs)),
                      "finufft_evaluations": sorted(set(sum(s["evaluations"] for s in p["finufft"]["passes"]) for p in pairs))})
    assert summary["all_acceptance_checks_pass"]
    result = {"passed": True, "endpoints_audited": len(audits), "direct_reference": audits,
              "max_audited_energy": max(r["energy"] for r in audits),
              "max_audited_S": max(r["max_S"] for r in audits),
              "max_audited_force_component": max(r["max_force_component"] for r in audits),
              "max_energy_difference_from_cpp": max(abs(r["energy_difference_from_cpp"]) for r in audits),
              "max_force_difference_from_cpp": max(abs(r["max_force_difference_from_cpp"]) for r in audits),
              "per_seed": table}
    (root/"independent_audit.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    with (root/"paired_runs.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(summary["optimization_pairs"][0]))
        writer.writeheader()
        writer.writerows(summary["optimization_pairs"])
    med_speed = summary["wall_speedup_median"]
    micro = summary["microbenchmarks"]
    lines = ["# Direct versus FINUFFT Type 1: optimization comparison", "",
             "All 30 optimizations passed: five seeds, three repetitions, two backends.", "",
             "## Conditions", "",
             "N=500, 2D square periodic box, number density=1, side=sqrt(500). "
             "299 independent half-plane mode pairs, chi=299/998=0.299599. "
             "The outer shell is partially selected to preserve the unit-test mode list; this is not a full-shell isotropic cutoff.", "",
             "Both backends reread the same saved initial.tsv for every run. Unit mode weights, zero energy shift, "
             "no repulsion. FINUFFT tolerance=1e-12. Existing RelaxStructure_NLOPT and objective code were compiled "
             "unchanged: three L-BFGS passes (scales 10000,1,0.01), memory100, max2000 evaluations/pass, "
             "original stopping defaults. All 90 passes returned NLopt success (1). "
             "Objective evaluations are reported; NLopt iterations are not exposed by this wrapper.", "",
             f"Two threads and two distinct physical cores (CPUs {settings['cpu_affinity']}); "
             "OMP_WAIT_POLICY=PASSIVE, OMP_PROC_BIND=false for both backends. Runs execute serially, "
             "alternating which backend runs first. Fresh sources are compiled at -O2; all timing includes "
             "the same lightweight counters. Standalone configuration validation and file output are outside optimizer wall time.", "",
             "Ground-state acceptance: direct total energy <=1e-11, max constrained S(k)<=1e-11, "
             "max absolute force component<=1e-7, no evaluation-limit stops, and no FINUFFT fallback.", "",
             "## Results", "",
             f"Median paired wall-time speedup: **{med_speed:.3f}x** "
             f"(range {summary['wall_speedup_range'][0]:.3f}-{summary['wall_speedup_range'][1]:.3f}x), "
             f"equivalent to {100*(1-1/med_speed):.1f}% less optimization time.", "",
             "| Seed | Direct median s | FINUFFT median s | Direct E (max) | FINUFFT E, direct reference (max) | Aligned RMS (max) |",
             "|---|---:|---:|---:|---:|---:|"]
    for row in table:
        lines.append(f"| {row['seed']} | {row['direct_seconds']:.4f} | {row['finufft_seconds']:.4f} | "
                     f"{row['direct_energy_max']:.3e} | {row['finufft_energy_max']:.3e} | {row['aligned_rms_max']:.3e} |")
    lines += ["", "Positions retain particle labels. Alignment removes periodic-image choices and the best "
              "global translation, but does not permute particles. Cartesian RMS difference range: "
              f"{summary['aligned_position_rms_range'][0]:.3e}-{summary['aligned_position_rms_range'][1]:.3e}. "
              "For these seeds the two methods reached numerically indistinguishable configurations; "
              "this is not a general uniqueness claim for stealthy ground states.", "",
              "### Evaluation microbenchmarks", "",
              "Each of the five initial and ten final configurations was evaluated by both backends, "
              "with five warm-up rounds and 30 interleaved timing samples. Table entries are medians "
              "across these 15 same-configuration medians, not a separately optimized input per backend.", "",
              "| Operation | Direct ms | FINUFFT ms | Median paired speedup |", "|---|---:|---:|---:|"]
    for field, title in (("rho_seconds","Density"), ("cached_force_seconds","Forces, cached density"), ("energy_force_seconds","Energy + forces")):
        dm=statistics.median(r[field]["direct_median"] for r in micro)*1000
        fm=statistics.median(r[field]["finufft_median"] for r in micro)*1000
        speed=statistics.median(r[field]["speedup"] for r in micro)
        lines.append(f"| {title} | {dm:.4f} | {fm:.4f} | {speed:.3f}x |")
    direct_evals=[sum(s["evaluations"] for s in p["direct"]["passes"]) for p in all_pairs]
    fast_evals=[sum(s["evaluations"] for s in p["finufft"]["passes"]) for p in all_pairs]
    lines += ["", f"Objective evaluations: direct {min(direct_evals)}-{max(direct_evals)}, "
              f"FINUFFT {min(fast_evals)}-{max(fast_evals)}. Force evaluation remains direct and limits overall speedup.", "",
              "## Independent saved-coordinate audit", "",
              "audit_comparison.py reads every saved endpoint and mode list and independently computes scalar "
              "direct sums using Python math.fsum; it calls neither C++ nor FINUFFT. "
              "It also recomputes periodic/aligned position differences from saved files.", "",
              f"All {len(audits)} endpoints passed. Maximum independent energy: {result['max_audited_energy']:.3e}; "
              f"maximum S(k): {result['max_audited_S']:.3e}; maximum force component: {result['max_audited_force_component']:.3e}. "
              f"Largest differences from C++ validation: energy {result['max_energy_difference_from_cpp']:.3e}, "
              f"max-force magnitude {result['max_force_difference_from_cpp']:.3e}.", "",
              "## Performance caveat", "",
              "The single-seed exploratory pilots are retained beside full/. The first accidentally selected "
              "SMT siblings, so it is not the primary benchmark. On distinct physical cores, the default "
              "OpenMP waiting behavior gave approximately 0.77x total speedup unbound and 0.99x with close binding; "
              "passive waiting gave 1.46x. This demonstrates runtime-setting sensitivity, not a universal "
              "FINUFFT speedup. Idle OpenMP workers competing with the DUCC0 thread pool is a plausible explanation; "
              "it was not established by a scheduler trace. Both backends use identical settings within each experiment.", "",
              "## Artifacts and scope", "",
              "settings.json records environment, affinity, binary/library hashes and thresholds. Per-seed JSON "
              "contains raw timing samples, pass statuses, counters and endpoint checks. Seed directories contain "
              "initial.tsv, modes.tsv, both backend endpoint TSVs and direct-reference structure factors. "
              "paired_runs.csv and independent_audit.json provide machine-readable summaries. "
              "The binary and exact build commands/source hashes are in ../pilot2/build/.", "",
              "This validates the base ShiftedCCPotential Type1 path with the existing CCOptimization optimizer. "
              "No application CLI switch, Type2 acceleration, five-component run, or feasibility sweep was added.", ""]
    (root/"REPORT.md").write_text("\n".join(lines))
    print(json.dumps({k:v for k,v in result.items() if k not in ("direct_reference","per_seed")}, indent=2))


if __name__ == "__main__":
    main()
