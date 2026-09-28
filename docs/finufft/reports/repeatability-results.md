# FINUFFT repeatability fix — completed September 27, 2026

Same-backend repeatability is restored in the tested five-component configurations. This is finite evidence for the tested machine, build, and settings, not a cross-platform guarantee or a claim that different backends select identical ground states.

## Production change

`CollectiveCoordinatePotential.h` now defaults FINUFFT to one internal thread, independently of `ParallelNumThread`. Other OpenMP work retained four threads in this experiment. The existing fixed-order contact-energy sum is preserved.

`SetFinufftNumThreads(int)` and `GetFinufftNumThreads()` expose the setting. The setter rejects values below one and invalidates cached collective coordinates. Explicit multithreaded overrides remain available without a repeatability guarantee. Backend selection, physical formulas, forces, and FINUFFT tolerance are unchanged; direct remains the default backend.

## Executed regression

- Five equal components, total N=2000, two dimensions, component chi approximately 0.1, total packing fraction 0.2.
- Seeds 20260928, 20260929, and 20260930; two runs per backend per seed: 12 optimizations.
- Initial positions and constrained modes were byte-identical across each matched set.
- Repeated endpoint files and optimization pass records were byte-identical within each backend for every seed.
- All 12 endpoints passed independent ground-state acceptance: energy <=1e-11, maximum structure factor <=1e-10, maximum force <=1e-7, and overlap <=1e-6.
- Fresh regression passed 101 required checks. Enabled and disabled FINUFFT API tests passed 75 and 50 checks respectively.
- Fixed-coordinate probe had zero energy and force range across 100 repetitions per backend. Its observation-only wrapper recorded 565 successful FINUFFT calls, all with nthreads=1; it did not modify call options.

### Optimization time in seconds

| Seed | Vanilla repeat 1 | Vanilla repeat 2 | FINUFFT repeat 1 | FINUFFT repeat 2 |
| --- | ---: | ---: | ---: | ---: |
| 20260928 | 6.665 | 7.436 | 6.015 | 5.978 |
| 20260929 | 6.674 | 6.625 | 6.149 | 5.892 |
| 20260930 | 6.084 | 6.061 | 5.023 | 5.033 |

Median optimization times: vanilla 6.645 s, FINUFFT 5.935 s. These are end-to-end optimization measurements with potentially different evaluation counts, not isolated transform benchmarks. Process timings and individual receipts are also retained.

Cross-backend endpoint differences remain. Removing run-to-run variability does not make the direct and approximate FINUFFT numerical trajectories identical.

## Independent audit

The requested multiagent skill separated Claude implementation from independent Codex critique and verification. Cycle 1 passed both gates; no second or third cycle was necessary. The verifier independently reran API tests and the fixed-coordinate probe and inspected all 12 optimization receipts, hashes, matching initial states, and acceptance criteria. It did not rerun the full optimization suite.

Two low-priority, nonblocking audit findings remain: the harness should check the return code of `ldd`, and its README should describe cached-rho assertions as exact value equality rather than bitwise equality. Actual dependency inspection succeeded. Unchanged forced failure/fallback paths were not newly exercised.

## Evidence and handoff

- [Full regression report](evidence-cycle1/REPORT.md)
- [Machine-readable results](evidence-cycle1/summary.json)
- [Build, source, library, and binary manifest](evidence-cycle1/manifest.json)
- [Independent code audit](audit-cycle1.md)
- [Independent verification](verification-cycle1.md)

Changes are uncommitted. The broader scientific sweep was not resumed; final service inspection confirmed `finufft-five-component-20260928.service` inactive with MainPID=0. Whitespace validation passed with the existing CRLF convention recognized. Unrelated worktree changes were preserved.
