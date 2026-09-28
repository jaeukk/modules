# Direct versus FINUFFT Type 1: optimization comparison

All 30 optimizations passed: five seeds, three repetitions, two backends.

## Conditions

N=500, 2D square periodic box, number density=1, side=sqrt(500). 299 independent half-plane mode pairs, chi=299/998=0.299599. The outer shell is partially selected to preserve the unit-test mode list; this is not a full-shell isotropic cutoff.

Both backends reread the same saved initial.tsv for every run. Unit mode weights, zero energy shift, no repulsion. FINUFFT tolerance=1e-12. Existing RelaxStructure_NLOPT and objective code were compiled unchanged: three L-BFGS passes (scales 10000,1,0.01), memory100, max2000 evaluations/pass, original stopping defaults. All 90 passes returned NLopt success (1). Objective evaluations are reported; NLopt iterations are not exposed by this wrapper.

Two threads and two distinct physical cores (CPUs [0, 2]); OMP_WAIT_POLICY=PASSIVE, OMP_PROC_BIND=false for both backends. Runs execute serially, alternating which backend runs first. Fresh sources are compiled at -O2; all timing includes the same lightweight counters. Standalone configuration validation and file output are outside optimizer wall time.

Ground-state acceptance: direct total energy <=1e-11, max constrained S(k)<=1e-11, max absolute force component<=1e-7, no evaluation-limit stops, and no FINUFFT fallback.

## Results

Median paired wall-time speedup: **1.450x** (range 1.402-1.506x), equivalent to 31.0% less optimization time.

| Seed | Direct median s | FINUFFT median s | Direct E (max) | FINUFFT E, direct reference (max) | Aligned RMS (max) |
|---|---:|---:|---:|---:|---:|
| 20260927 | 1.0783 | 0.7398 | 5.113e-22 | 7.817e-22 | 1.154e-11 |
| 20260928 | 1.0443 | 0.7069 | 4.400e-22 | 1.005e-21 | 7.681e-12 |
| 20260929 | 1.0469 | 0.7375 | 2.885e-22 | 7.099e-22 | 3.438e-12 |
| 20260930 | 1.0468 | 0.7355 | 9.260e-22 | 9.975e-22 | 8.094e-12 |
| 20261001 | 1.0883 | 0.7469 | 1.356e-21 | 1.620e-21 | 3.275e-11 |

Positions retain particle labels. Alignment removes periodic-image choices and the best global translation, but does not permute particles. Cartesian RMS difference range: 3.318e-12-3.275e-11. For these seeds the two methods reached numerically indistinguishable configurations; this is not a general uniqueness claim for stealthy ground states.

### Evaluation microbenchmarks

Each of the five initial and ten final configurations was evaluated by both backends, with five warm-up rounds and 30 interleaved timing samples. Table entries are medians across these 15 same-configuration medians, not a separately optimized input per backend.

| Operation | Direct ms | FINUFFT ms | Median paired speedup |
|---|---:|---:|---:|
| Density | 1.8017 | 0.5622 | 3.191x |
| Forces, cached density | 2.0397 | 2.2769 | 0.917x |
| Energy + forces | 4.0535 | 2.7873 | 1.451x |

Objective evaluations: direct 245-253, FINUFFT 241-258. Force evaluation remains direct and limits overall speedup.

## Independent saved-coordinate audit

audit_comparison.py reads every saved endpoint and mode list and independently computes scalar direct sums using Python math.fsum; it calls neither C++ nor FINUFFT. It also recomputes periodic/aligned position differences from saved files.

All 30 endpoints passed. Maximum independent energy: 1.621e-21; maximum S(k): 5.052e-22; maximum force component: 7.824e-12. Largest differences from C++ validation: energy 1.919e-25, max-force magnitude 1.834e-14.

## Performance caveat

The single-seed exploratory pilots are retained beside full/. The first accidentally selected SMT siblings, so it is not the primary benchmark. On distinct physical cores, the default OpenMP waiting behavior gave approximately 0.77x total speedup unbound and 0.99x with close binding; passive waiting gave 1.46x. This demonstrates runtime-setting sensitivity, not a universal FINUFFT speedup. Idle OpenMP workers competing with the DUCC0 thread pool is a plausible explanation; it was not established by a scheduler trace. Both backends use identical settings within each experiment.

## Artifacts and scope

settings.json records environment, affinity, binary/library hashes and thresholds. Per-seed JSON contains raw timing samples, pass statuses, counters and endpoint checks. Seed directories contain initial.tsv, modes.tsv, both backend endpoint TSVs and direct-reference structure factors. paired_runs.csv and independent_audit.json provide machine-readable summaries. The binary and exact build commands/source hashes are in ../pilot2/build/.

This validates the base ShiftedCCPotential Type1 path with the existing CCOptimization optimizer. No application CLI switch, Type2 acceleration, five-component run, or feasibility sweep was added.
