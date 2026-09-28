# Vanilla versus FINUFFT endpoint agreement: unresolved

September 27, 2026. This task addresses different backends starting from the same initial state, not the already established repeatability within either backend.

## Outcome

No tested correction restored positional agreement. Production code was not changed in this task. The diagnostic evidence strongly supports numerical trajectory sensitivity in the current three-pass L-BFGS workflow: tiny evaluator differences remain tiny when measured at identical coordinates, while the separately optimized endpoints differ substantially.

This is a bounded conclusion, not proof that agreement is impossible, that every integration bug has been excluded, or that the dynamics are chaotic. Changing the optimizer is a separate decision: it can change the ground state selected by vanilla as well as by FINUFFT. An optional alternative-optimizer experiment was proposed to the user; no alternative was implemented or run pending that choice.

## Conditions and checks

- N=2000, two dimensions, five equal components, component chi=80/798, total packing fraction=0.2, identical radii.
- Four OpenMP threads outside FINUFFT, one internal FINUFFT thread, fixed-order contact-energy summation retained.
- Seeds 20260928, 20260929, 20260930. Initial coordinates and mode files were matched by seed. Baseline traced endpoints reproduce the preceding verified benchmark.
- Seventeen optimizations: six baseline runs, three vanilla-only perturbation controls, two tighter-FINUFFT-tolerance runs, and six shorter-L-BFGS-history runs.
- All resulting states passed the independent energy, structure-factor, force, and overlap criteria. History-1 FINUFFT hit the 2000-evaluation limit in its first pass (status 5), then completed the later passes (status 1); it was not uninterrupted convergence.
- Endpoint agreement has a separate provisional engineering gate: maximum fixed-label, periodic, globally translation-aligned fractional displacement below 1e-7. Every tested paired backend comparison failed this gate. This tolerance was not specified by the user.

## Paired results

Displacements below are mean Cartesian particle distances, not fractional-coordinate distances. Times are optimization wall times with trace-writing overhead, not clean performance benchmarks. The direct timing in tighter-epsilon rows is the same baseline run, not an additional run.

| Seed / control | Vanilla seconds | FINUFFT seconds | Mean endpoint displacement |
| --- | ---: | ---: | ---: |
| 20260928, default history 100, eps=1e-12 | 6.580 | 5.914 | 0.430886 |
| 20260929, default history 100, eps=1e-12 | 6.637 | 5.918 | 0.410684 |
| 20260930, default history 100, eps=1e-12 | 5.972 | 5.064 | 0.286372 |
| 20260928, eps=1e-14 | 6.580 | 5.772 | 0.423731 |
| 20260928, eps=1e-15 | 6.580 | 6.091 | 0.413784 |
| 20260928, history=1 | 10.150 | 12.568 | 0.182479 |
| 20260928, history=5 | 6.791 | 6.704 | 0.217490 |
| 20260928, history=10 | 6.294 | 5.647 | 0.246795 |

Tighter epsilon was applied by a diagnostic wrapper to the actual FINUFFT call, bypassing the production eligibility guard while leaving the accepted integer mode mapping unchanged. All recorded calls succeeded, used one internal thread, and the optimization backend counters showed actual FINUFFT use. These are diagnostic controls, not a proposed production guard change or a direct-sum fallback disguised as FINUFFT.

## Trajectory evidence

For seed 20260928 at the identical first input, the optimizer-space gradients differed by at most 1.307e-16 at epsilon=1e-12. The recorded fractional callback-coordinate difference exceeded 1e-8 at callback 57 and 0.01 at callback 149.

These are callback indices, **not accepted-iteration numbers**. Different line-search decisions can cause equal-index callbacks to represent different stages; no iteration-by-iteration amplification rate has been established.

A diagnostic that repeatedly added 1e-17 to the first gradient component in vanilla, without changing its returned energy, produced mean endpoint shifts of 0.392384, 0.372276, and 0.284647 across the three seeds. This deliberately inconsistent-gradient control demonstrates sensitivity without invoking FINUFFT. It is neither a physical correction nor a faithful model of FINUFFT error. Perturbed-run optimization times were 6.076, 7.250, and 5.929 seconds.

Removing independent translations of each species barely changed the original discrepancy, ruling out simple component drift as an explanation for these samples.

## Frozen-coordinate replay

Both backends evaluated exactly the same saved configurations from both trajectories: 66 first-pass samples plus 24 later-pass samples, covering all three seeds and all three optimization passes. This replay performs no optimization and uses fresh potential objects, direct Fourier/all-pair reference calculations, and actual-backend/status checks for all five species.

Across all 90 snapshots:

| Quantity | Maximum difference |
| --- | ---: |
| Vanilla versus FINUFFT total energy | 4.2633e-13 |
| Physical force component | 3.4564e-14 |
| Complex collective coordinate | 7.0089e-11 |
| Matching backend replay versus saved energy | 0 |
| Matching backend replay versus saved optimizer gradient | 0 |

All replay gates passed with successful FINUFFT use. Zero saved-value errors establish exact numeric equality, not an asserted bytewise representation comparison. The evidence finds no state-dependent mapping/force discrepancy at the sampled configurations, despite the differing optimizer endpoints.

## Audits and remaining decision

The multiagent skill separated Claude-authored diagnostic code from independent Codex criticism and verification. The initial oversized implementation attempt was stopped after more than 12 minutes without producing files; it left no code changes. Narrower managed writes produced the trace and replay tools.

- [Audit cycle 1](audit-cycle1.md): diagnostic PASS; endpoint objective FAIL. Requested later same-coordinate replay.
- [Audit cycle 2](audit-cycle2.md): replay method and bounded conclusion PASS; endpoint objective still FAIL. Reviewed the first 66 replay samples; the subsequent 24 samples used the same unchanged replay code.
- [Independent executed verification](verification.md): reran early- and later-pass replay selections, checked repeat equality, all 17 optimization receipts, all 90 saved replay results, matched initial states, baseline endpoint reproduction, and 26 source/object/binary hashes. Diagnostic evidence PASS; endpoint objective FAIL.
- [Matched-control pairing clarification](verification-pairing.md): independently recomputed the five epsilon/history comparisons with the correct corresponding direct run. This supersedes the first verification note's mixed-control aggregate of 0.06151. The maximum fractional norm among these five correctly matched controls is 0.0429480; all fail 1e-7. The follow-up script checked IDs and species (its column-1 assertion was mislabeled as radii); the conductor additionally checked IDs, species, and actual radius column 4 for all five pairs and found exact matches.

Two of the allowed three audit cycles have been used. The third remains available for a candidate correction after the optimizer-choice decision. No failed endpoint gate has been converted into success by accepting low energy alone, copying endpoints, quantizing values, or replacing FINUFFT with direct sums.

## Reproducible artifacts

- [Trace and control analysis](evidence1/analysis.json), [additional controls](evidence1/controls.json)
- [First-pass replay summary](evidence1/replay_summary.json), [later-pass replay summary](evidence1/replay_later_summary.json)
- [Build and source/object hashes](evidence1/manifest.json), [replay build hashes](evidence1/replay_manifest.json)
- Per-run `*.receipt.json`, coordinate/mode files, and native-double `trace.bin` files are retained in `evidence1/`. Each optimization receipt records the exact command, return code, thread/tolerance observations, elapsed optimization time, and validation results.
- Replay receipts ending in `.replay2.receipt.json` separately capture stdout and stderr. `direct.replay.receipt.json` is superseded: its combined streams interleaved and caused a parsing failure, so the read-only replay was rerun with separate capture.

The trace builds reused previously verified core objects after checking their source hashes; only diagnostic entry points were newly compiled. Shared-object hashes and commands are recorded. Existing production edits were preserved, all new work remains uncommitted, and the broader scientific sweep remains stopped.
