**No blocking findings.** The reviewed fix resolves the round-1 slowdown in the tested cases.

**FOLLOW-UP — Validate the crossover on Nurion.**

- **Location:** `potential/CollectiveCoordinatePotential.h:399–429`.
- **Trigger:** FINUFFT selected on Xeon Phi 7250 KNL or Xeon 6148 Skylake, with a workload near the gate’s acceptance threshold.
- **Impact:** Machine-specific constants could accept a transform that is slower than the baseline. The 1.3 margin does not guarantee portability.
- **Judgment:** FOLLOW-UP, not blocking. A slower CPU alone does not establish a regression: both paths change cost, and assuming ideal direct-loop parallel scaling makes the gate conservative as thread count increases. I have neither Nurion measurements nor a defensible quantitative case showing a beyond-noise slowdown there. The documented production-node calibration remains appropriate.

The implementation supports the determinism claim: gate decisions use configuration, constraints, and settings, with no timing feedback. Declined species retain the original fused `Force(i)` arithmetic. Splitting mode mapping from angle filling preserves the inspected Type 1 calculations.

I built a disposable harness using `bench_finufft.cpp`’s packing and potential setup. It compared current `AllForce` with the reconstructed prior parallel `Force(i)` loop, including contact forces and using the same cached FINUFFT rho. Execution order alternated; measurements discarded 10 warmups and used 100 samples.

Results with default OpenMP waiting:

| Case | Type 2 species | Current / prior median time | Bitwise equal to prior |
|---|---:|---:|---|
| 3D, N=64, S=1, χ=0.3, T=8 | 0 | 1.010× | Yes |
| 3D, N=320, S=10, χ=0.3, T=8 | 0 | 1.004× | Yes |
| 3D, N=1000, S=10, χ=0.3, T=8 | 0 | 1.019× | Yes |
| 2D, N=2000, S=5, χ=0.4, T=4 | 5 | 0.384× | No, as expected |
| 2D, N=2000, S=5, χ=0.4, T=8 | 5 | 0.500× | No, as expected |

Every case repeated bitwise across all 110 evaluations. A passive-wait run also showed no declined-path regression: ratios were 0.998×, 0.985×, and 0.921×; accepted cases were 0.391× and 0.399×. Timing variation limits precision, but the original 3–12× slowdown is eliminated and the claimed accepted-case speedups are supported.

Commands and results:

| Command | Result |
|---|---|
| `make -C apps/CCOptimization test-finufft-forces` | Passed |
| `make -C apps/CCOptimization test-finufft` | Passed |
| `bash potential/tests/finufft_type1/run_tests.sh` | All variants passed: 40,977 / 40,999 / 40,999 checks |
| Forced recompilation of both app test sources and `MultiStealthyPotential.cpp`, then reran both executables | Both passed |
| `make -C apps/CCOptimization -f review.mk USE_FINUFFT=1 review_gate.out` | Harness built successfully |
| `OMP_DYNAMIC=false OMP_WAIT_POLICY=PASSIVE apps/CCOptimization/review_gate.out` | Declines, repeatability, and timings above |
| `OMP_DYNAMIC=false apps/CCOptimization/review_gate.out` | Default-wait results above |

**Review limitation:** This copy has no Git metadata. Scope was restricted to the supplied fix diff and its affected paths; baseline comparisons reconstruct the prior loop rather than execute a checkout of `c352d19`. Nurion performance was not measured.