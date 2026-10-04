# Paired optimization comparison

Run from the cloned `modules/` directory. Choose a new output directory; existing experiments are never overwritten.

```sh
python3 potential/tests/finufft_type1/run_comparison.py --output /absolute/path/to/new-results
python3 potential/tests/finufft_type1/audit_comparison.py /absolute/path/to/new-results
```

The runner compiles fresh dependency objects into the experiment directory (two compiler jobs), then runs experiments serially. FINUFFT (`FINUFFT_PREFIX`), GSL, and NLopt default to `/home/jaeukk/install`. It does not rebuild or modify application objects. `--binary PATH` reuses an existing benchmark executable, whose hash is recorded.

Default design: N=500, 2D, number density1, 299 independent mode pairs, chi=299/998, five fixed seeds, three repetitions per backend, two threads, two distinct physical cores, FINUFFT tolerance1e-12. Both runs reload one saved starting configuration. Backend order alternates. The same unmodified `RelaxStructure_NLOPT` used by CCOptimization performs three L-BFGS passes. A link-time wrapper records return codes and objective counts without changing the algorithm. Each pass allows2000 objective evaluations; original stopping settings are preserved.

Acceptance requires direct-reference energy<=1e-11, maximum constrained structure factor<=1e-11, maximum force component<=1e-7, successful optimizer return codes without budget exhaustion, and no FINUFFT fallback. Ground-state criteria were set before the pilot. Endpoint positions are compared with fixed particle labels, periodic boundaries, and a least-squares global translation. No permutation matching or claim of ground-state uniqueness is made.

`audit_comparison.py` independently reads every saved endpoint and recomputes direct densities, energy, and forces with Python scalar arithmetic and `math.fsum`, without linking C++ or FINUFFT. It checks reported displacements, saves `independent_audit.json` and `paired_runs.csv`, and generates `REPORT.md` for the default full experiment.

Thirty timing samples after five warm-ups compare density, cached-density forces, and energy-plus-forces on each initial configuration and both first-repetition endpoints. Optimization timings exclude warm-up, endpoint validation and output; density and force counters are applied equally to both methods. Reported counts are objective evaluations, not NLopt iterations.

## Threading caveat

The benchmark defaults to `OMP_WAIT_POLICY=PASSIVE` and `OMP_PROC_BIND=false`, equally for both backends, with process affinity restricted to two distinct physical cores. The installed FINUFFT uses DUCC0 and OpenMP; timing was sensitive to idle-wait policy. Single-seed pilots with default OpenMP waiting were slower than the passive-waiting runs and remain in the experiment archive. The reported speedup is conditional on these settings, not a universal performance claim.

`--cpus CPU CPU`, `--omp-proc-bind close`, and `--omp-wait-policy ACTIVE` are explicit diagnostic overrides. Neighboring CPU numbers can be SMT siblings; the default selects distinct physical core IDs.

## Recorded experiment

The completed study is stored at:

`_multiagent/tasks/2026-09-27-finufft-comparison/full/` relative to the outer `multi-stealthy` project.

All30 endpoints passed an independent direct-sum audit. Median paired total speedup was1.450x; median density speedup was3.191x. Aligned position RMS differences ranged from3.318e-12 to3.275e-11. See the generated report and raw data for settings, per-seed measurements, and limitations.

This exercise does not add a production CLI switch, Type2 force acceleration, or a five-component feasibility study.
