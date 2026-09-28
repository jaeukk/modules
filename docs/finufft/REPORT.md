# FINUFFT integration: consolidated findings and release decision

September 28, 2026. This document supersedes earlier reports for current implementation status. Historical report snapshots are preserved in [reports/manifest.json](reports/manifest.json), including source paths and hashes. Their relative links refer to original evidence folders and may not resolve within these snapshots. Dates, "running," "uncommitted," and incomplete-audit statements in old reports describe their original snapshots, not current status.

## Release decision

FINUFFT exists to accelerate collective-coordinate evaluation. The approximately 90-fold slower fixed-parameter optimizer is **rejected**, even though it produced matching five-component ground states in one paired diagnostic. None of the factor investigation's alternative optimizers, force-sum variants, tolerance changes, or contact-disabled controls is adopted.

Production keeps the optional FINUFFT Type 1 density path and the existing `RelaxStructure_NLOPT()` three-pass L-BFGS procedure. Force evaluation remains direct. Direct summation remains the default backend. `StructureOptimization.cpp` and `.h` were not modified by the factor investigation, so there is no production optimizer patch to revert.

The remaining limitation is explicit: with five components and contact repulsion, matching original initial states can produce distinct ground-state configurations across Direct and FINUFFT, although both pass independent physical checks. Same-backend repeatability is a separate property and is restored in the tested settings. This release does not claim cross-backend configuration identity.

## What is retained and what is rejected

| Change or control | Disposition | Demonstrated effect / limitation |
| --- | --- | --- |
| Optional FINUFFT Type 1 `GetRho()` implementation | Retain | Accelerates densities on eligible reciprocal grids; reports actual backend and fallback status. No Type 2 force acceleration. |
| Species bookkeeping, cache initialization, radius/cutoff correctness fixes | Retain, already committed | Required correctness fixes, not alternative optimization physics. |
| Fixed-order contact-energy sum | Retain | Restored vanilla same-backend repetition in tested five-component runs. Parallel per-particle work followed by particle-index summation. |
| Independent FINUFFT thread setting, default 1 | Retain | Restored FINUFFT same-backend repetition in tested runs. Other OpenMP work still uses four threads; explicit multithreaded FINUFFT remains available without that guarantee. |
| Parameterized regression harnesses and thread observer | Retain | Verify actual backend use, optional/disabled builds, matching starts, endpoint checks and repetitions. |
| Changed L-BFGS scales/history, tighter transform accuracy | Diagnostic only | Did not restore cross-backend endpoint agreement. Production defaults unchanged. |
| Sorted, long-double contact-force summation | Retired | No cross-backend fix; was never production code. Do not confuse this with the retained contact-energy sum. |
| Contact-disabled five-component control | Diagnostic only | Closely matched endpoints, but changes the objective and is not a full-model ground state. |
| Truncated Newton / adaptive spectral descent | Retired | Failed agreement and/or convergence within their budgets. |
| Fixed-parameter momentum optimizer | Rejected and retired | Matching endpoints for one seed, but approximately 550/501 seconds versus seconds for L-BFGS. |

All active factor-only source files and runner entry points were removed from the research task folder; recoverable source and binaries are archived in `_multiagent/tasks/finufft-factors/retired/`. The measured evidence remains. These experimental sources are not included in this code release. Previously existing correctness fixes were not reverted.

## Evidence in chronological order

### One component: agreement and acceleration

N=500, 2D, chi=299/998, five seeds and three repetitions per backend, no contact: all 30 optimizations passed. Median paired speedup was 1.450x. Direct median times ranged 1.0443–1.0883 seconds; FINUFFT 0.7069–0.7469 seconds. Aligned Cartesian RMS differences ranged 3.318e-12–3.275e-11. Density speedup was 3.191x; unchanged direct force work limited overall speedup. These results depend on the tested affinity and passive OpenMP waiting settings, not a universal performance claim. [Full benchmark](reports/one-component-benchmark.md).

**Correction to an overly broad explanation:** one-component cases *with contact repulsion* were also tested, using three matched seeds at two and four threads. Mean cross-backend differences stayed 1.53e-11–5.78e-11. Contact alone therefore does not explain why the five-component sample is sensitive. Component count, chi, relative stiffness and trajectory effects have not been individually isolated. [Contact study](reports/contact-results.md).

### Five components: two distinct reproducibility issues

N=2000, 400 particles per component, equal radii, chi per component=80/798, total packing fraction=0.2. Initial state and mode matching were checked, not assumed from seed labels.

1. Original parallel contact-energy reduction caused vanilla run-to-run differences. Fixed-order energy summation removed them in the tested runs.
2. FINUFFT multithreaded summation still caused run-to-run variation. Setting FINUFFT's internal thread count to one independently of outer OpenMP restored tested same-backend repetition. This does not make its arithmetic identical to direct summation.
3. Cross-backend endpoints still differ, with mean Cartesian displacements approximately 0.430886, 0.410684 and 0.286372 for seeds 20260928–20260930. Both endpoints pass physical checks. [Thread-control study](reports/single-thread-control.md), [production repeatability fix](reports/repeatability-results.md).

The original stopped sweep is historical and not restarted. Its early endpoint differences cannot be attributed solely to FINUFFT because same-backend variability was present then. [Original snapshot](reports/five-component-original.md), [pairing audit](reports/five-component-pairing.md).

### Same-coordinate and trajectory checks

Across 90 saved-coordinate replays, maximum Direct/FINUFFT differences were 4.2633e-13 in total energy, 3.4564e-14 in force components, and 7.0089e-11 in complex collective coordinates. Matching-backend saved energies and gradients were reproduced exactly. This supports small evaluator differences followed by optimization sensitivity; it does not prove every possible bug excluded or identify one unique line-search decision. [Cross-backend audit](reports/cross-backend-results.md).

For seed 20260928, all 2291 baseline callbacks were replayed with separate CC/contact energies. Divergence arose well before ground-state energies. Looser stopping targets did not resolve it; removing individual component drifts reduced mean displacement only about 0.067%. These are callback counts, not accepted L-BFGS iterations. [Energy/drift/tolerance study](reports/energy-diagnosis.md).

The endpoint pair also failed a necessary geometric-identity condition: complete periodic pair-distance spectra differed, even ignoring component labels. Translation, lattice-compatible rotation/inversion and particle relabeling cannot explain the discrepancy. This invariant is a rejection test, not a sufficient identity test when spectra agree. [Zhang-based identity assessment](reports/geometric-identity.md).

### Four factor-audit cycles: evidence, not a release fix

Coordinate scaling, fixed-order higher-precision contact forces, truncated Newton and adaptive spectral descent did not resolve the full-model comparison. Without contact, five-component endpoints matched to mean 1.78e-12, but the states were not ground states of the original full model.

Fixed-parameter momentum gave mean difference 3.80e-12 and maximum fractional difference 4.03e-12; both independent final energies were approximately 1.0e-13. This used each backend independently from the same original seed and 124043 evaluations per backend, including a checkpoint restart with velocity reset. It was tested for one original seed; a complete end-to-end repetition of this alternative was not performed.

L-BFGS refinement amplified discrepancies again. Even starting from exactly the same late checkpoint gave mean difference 7.78e-5. Repeated refinement was exact within each backend. The findings demonstrate optimizer sensitivity, not a new threading race, and not an implementation bug proved in `RelaxStructure_NLOPT()`. The slow alternative is explicitly rejected. [Retired study](reports/retired-factors.md), [four-cycle ledger](reports/retired-factor-audits.md).

## Timing ledger: do not pool unlike experiments

| Study / settings | Vanilla seconds | FINUFFT seconds | Interpretation |
| --- | ---: | ---: | --- |
| One component with contact, fixed sum, two threads (median) | 1.313 | 0.917 | Matched ground states in tested seeds |
| One component with contact, fixed sum, four threads (median) | 0.718 | 0.678 | Matched ground states in tested seeds |
| Five components, fixed contact sum, FINUFFT four threads (median) | 7.142 | 10.234 | FINUFFT repeats still variable |
| Five components, retained one-thread FINUFFT production fix (median) | 6.645 | 5.935 | Same-backend repeatability restored; cross-backend states differ |
| Factor audit: single Cartesian-scale L-BFGS, seed 20260928 | 5.981 | 5.353 | Both ground states; mean difference 0.358722 |
| Rejected fixed-momentum two-stage workflow, same seed | 550.247 | 501.442 | Endpoint agreement at unacceptable cost |

These are optimization wall times, not physical simulation time. Different studies use different seeds, thread settings, optimizer paths and instrumentation; timings cannot be interpreted as a single controlled before/after speedup. Full reports preserve individual run tables.

## Release validation

A fresh build and the retained repeatability regression passed all 101 required checks, including 12 optimizations, independent endpoint checks and byte-identical within-backend repetitions. Release-run median times were 6.794 seconds vanilla and 6.102 seconds FINUFFT. Cross-backend mean differences remained 0.4309, 0.4107 and 0.2864 across the three seeds. See [release validation](RELEASE_VALIDATION.md). This validates the retained fast implementation, not the retired factor harness. No broader scientific sweep is authorized or started by this handoff.
