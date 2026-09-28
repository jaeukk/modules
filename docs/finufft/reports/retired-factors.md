# Four-cycle investigation of cross-backend configuration agreement

> RETIRED September 28, 2026: the user rejected the runtime tradeoff. No experimental optimizer or force-sum variant is adopted. Active diagnostic source/entry points were removed and preserved in `retired/diagnostic-source.tar.gz`; binaries were moved to `retired/`. Results below are historical evidence, not a current recommended workflow. The reproduction commands below require restoring the archive and are intentionally no longer active. Production `RelaxStructure_NLOPT()` was never modified by this investigation.

## Outcome

All four audit cycles are complete. A fixed-parameter momentum workflow passes both ground-state and cross-backend position criteria for the one tested original initial state. Mean Cartesian endpoint difference is 3.79989e-12; maximum fractional difference is 4.02850e-12, versus the unchanged 1e-7 gate. This is a demonstrated diagnostic alternative, not a production fix or a general guarantee across seeds, chi values, or system sizes. The original L-BFGS discrepancy remains.

## Scope and interpretation

One original initial state: seed 20260928, N=2000, 2D square cell of side sqrt(2000), five equal components, chi per component 80/798, total packing fraction 0.2. Four OpenMP threads and one internal FINUFFT thread; FINUFFT epsilon 1e-12. The existing fixed-order contact-energy sum remains enabled. Tests use genuine FINUFFT density evaluations, with call counters excluding silent fallback. Original initial coordinates are checked as arrays, not merely by seed labels.

The distinction is between deterministic repetition within one backend and agreement between two different floating-point evaluations of the same mathematical objective. The former already works; the latter is the issue tested here. No physical potential, ground-state acceptance criterion, or production source was changed. Experimental algorithms and force-order controls live only in this diagnostic folder.

Physical acceptance requires independently evaluated total energy <=1e-11, maximum constrained structure factor <=1e-10, maximum force component <=1e-7, and maximum relative overlap <=1e-6. The separate position gate is maximum periodic fixed-ID fractional Euclidean displacement <=1e-7. Mean distances below are Cartesian. This is stricter than searching for agreement only after translation or particle relabeling.

## Completed controls

| Cycle / control | Mean endpoint displacement | Both full-model ground states? | Position gate? |
| --- | ---: | --- | --- |
| 1: L-BFGS, coordinate scale 10000, first pass only | 0.430880 | No | Fail |
| 1: L-BFGS, scale 1 | 4.443408 | Yes | Fail |
| 1: L-BFGS, Cartesian scale | 0.358722 | Yes | Fail |
| 2: Same five-component system, contact disabled | 1.77986e-12 | No; intentionally different objective | Pass |
| 2: Full contact, sorted/long-double force accumulation | 0.362607 | Yes | Fail |
| 3: Preconditioned restarted truncated Newton | 0.745449 | No; evaluation budget reached | Fail |
| 3: Spectral steepest descent with Armijo search | 0.110246 | No; evaluation budget reached | Fail |
| 4: Fixed-parameter momentum, 40000 evaluations | 3.68697e-12 | No; intermediate state | Pass |
| 4: L-BFGS refinement of the preceding states | 6.69974e-5 | Yes | Fail |
| 4: L-BFGS from an exactly common late checkpoint | 7.78171e-5 | Yes | Fail |
| 4: Fixed-parameter continuation to ground accuracy | 3.79989e-12 | Yes | Pass |

The scale-10000 row reproduces only the old first pass, not the earlier full three-pass procedure. The no-contact control is not a proposed physical solution: it has contact energy about 134.19 when evaluated under the original full objective.

## What these controls establish

Coordinate scaling matters but does not remove the discrepancy. The five-component collective-coordinate problem without contact agrees extremely closely across backends. Sorting and increasing the precision of contact-force accumulation does not fix full-model agreement.

The clearest control is in cycle 4: fixed-parameter momentum preserves close endpoint agreement after 40000 evaluations, but L-BFGS refinement of those nearly identical states amplifies the displacement by roughly seven orders of magnitude. Repeating the refinement produces byte-identical initial/final files and energy histories within each backend. This is evidence of cross-backend sensitivity of the search process, not a newly observed same-backend race.

A separate late-stage control starts both backends from exactly the same vanilla momentum checkpoint. L-BFGS still yields mean displacement 7.78171e-5 and maximum fractional displacement 2.12367e-5. Thus the failed refinement cannot be attributed solely to its two slightly different incoming checkpoints. That common-checkpoint experiment is only a diagnostic; the successful fixed-parameter workflow uses each backend's own states throughout.

The initial state contains 786 overlapping pairs. The closest separation is 0.0105981 versus a contact diameter of 0.504627. The contact pair's radial curvature is +7.85398 and its most negative transverse curvature is -366.112 in relative Cartesian coordinates. These pair-level calculations show strong local nonconvexity; they are not full-system Hessian eigenvalues and do not identify a particular L-BFGS line-search branch as the unique cause. We have not proved a general stability theorem or uniqueness of any ground state.

## Fixed-parameter continuation: passed for this pair

The final part of cycle 4 continues each backend's own close checkpoint with velocity reset to zero in both runs. The update uses velocity retention 0.9 and Cartesian force coefficient 0.02, with no line search or curvature-history estimate. It is not physical molecular dynamics. The continuation allows at most 100000 evaluations and stops at native energy <=1e-13; all independent physical criteria and the position gate pass. The earlier 40000-step run used an unreached 1e-16 target. No endpoint is copied from the other backend and no direct-sum correction is applied to FINUFFT.

The continuation inputs differ slightly because they are each backend's own outputs. The two complete workflows nevertheless start from exactly the same original configuration. Saved `seed_initial.tsv`, `initial.tsv`, restart receipts/hashes, and runner assertions make this distinction explicit.

| Quantity | Vanilla | FINUFFT |
| --- | ---: | ---: |
| First-stage evaluations | 40000 | 40000 |
| Continuation evaluations | 84043 | 84043 |
| Total evaluations | 124043 | 124043 |
| First-stage seconds | 180.137 | 166.161 |
| Continuation seconds | 370.110 | 335.281 |
| Total optimization seconds | 550.247 | 501.442 |
| Independent total energy | 9.99800550836e-14 | 9.99800556612e-14 |
| Independent contact energy | 1.81071e-18 | 1.81071e-18 |
| Maximum force component | 9.20358e-10 | 9.20358e-10 |
| Maximum constrained structure factor | 1.39933e-12 | 1.39933e-12 |
| Maximum relative overlap | 2.17003e-10 | 2.16972e-10 |

FINUFFT was used for all 620215 density evaluations in its successful two-stage workflow. Each component's mean Cartesian endpoint discrepancy is between 2.47e-12 and 5.98e-12. No translation, rotation, inversion, or particle relabeling is needed to pass the position criterion.

Cost: approximately 9.17 minutes for vanilla and 8.36 minutes for FINUFFT, versus about 6 seconds for the Cartesian-scale L-BFGS ground-state test. The stable alternative is roughly 90 times slower in this diagnostic. Its speed has not been tuned, and the diagnostic timings are not clean benchmarks. Changing the default production optimizer would therefore be a material tradeoff, not an automatic bug fix.

### Reproduce the passing workflow

Run from the project root. Use fresh output names because the existing directories are protected against overwrite.

```sh
python _multiagent/tasks/finufft-factors/build.py
python _multiagent/tasks/finufft-factors/run.py NEW_stage1 3 44.721359549995796 1 40000 1
python _multiagent/tasks/finufft-factors/run.py NEW_finish 3 44.721359549995796 1 100000 1 NEW_stage1
```

The current diagnostic's fixed-descent stop target is 1e-13. This does not change the first stage here because its residual energy at 40000 evaluations is much larger. The SCALE and ENERGY_SCALE arguments apply to NLopt modes; the custom momentum method explicitly uses Cartesian coordinates and the fixed coefficients above.

## Validation and reproducibility

- Experimental ordered forces were checked at the original configuration and baseline callback 74. Maximum difference from original forces: 8.88178e-16. Maximum difference from the independent reference: 1.58673e-12. Maximum sampled central finite-difference error: 1.33834e-9.
- Builds verify hashes of reused production sources and object files. Later build manifests retain the complete diagnostic source text and binary hash. Receipts record the binary actually used; earlier build hashes are retained in receipts even where a full source-text snapshot was not saved.
- All runs retain initial/final positions and every objective evaluation's energy and gradient norm. Long diagnostic runs do not save all intermediate particle positions. The earlier dedicated trajectory investigation does retain its full particle histories.
- These are main-agent audits, not independent external-agent reviews. The complete fixed-parameter workflow has not been repeated end to end in this investigation; the L-BFGS refinement was repeated.
- Timings are diagnostic optimization wall times with logging and shared-machine conditions, not an isolated performance benchmark.

Files: [audit ledger](audit-cycles.md), [machine-readable analysis](analysis.json), [verification](verification.json), [probe source](probe.cpp), [build](build.py), [paired runner](run.py). Each `cycle*/summary.json` links its settings to retained command receipts and endpoint files. Output directories refuse overwrite.
