**BLOCKING — Unconditional Type 2 dispatch slows ordinary small-species workloads.**

- **Location:** `apps/CCOptimization/MultiStealthyPotential.cpp:480`; `potential/CollectiveCoordinatePotential.h:614`.
- **Trigger:** FINUFFT selected, 3D unit-density random packing, χ = 0.3, 8 OpenMP threads, N = 320 with 10 species.
- **Impact:** Every species performs a separate Type 2 transform regardless of workload size. Transform overhead makes `AllForce` substantially slower than the prior parallel per-particle force loop.
- **Evidence:** A disposable timing harness reused `bench_finufft.cpp`’s packing/potential setup and compared current `AllForce` against the prior parallel `Force(i)` loop, including contact forces, using the **same cached FINUFFT rho**. It alternated execution order, discarded 10 warmups, and reported medians over 100 samples.

Commands:
```sh
make -C apps/CCOptimization -f review.mk USE_FINUFFT=1 review_force_timing.out
OMP_DYNAMIC=false OMP_WAIT_POLICY=PASSIVE apps/CCOptimization/review_force_timing.out
```

Output after rebuilding:
```text
N=64 S=1 modes=61 threads=8 Type2_us=728.662 prior_loop_us=220.87 ratio=3.29905
N=320 S=10 modes=28 threads=8 Type2_us=8402 prior_loop_us=714.802 ratio=11.7543
N=1000 S=10 modes=89 threads=8 Type2_us=8959.81 prior_loop_us=3229.01 ratio=2.77479
```

A preceding passive-wait run showed slowdowns of 4.38×, 12.06×, and 3.69× respectively. Maximum force differences were below `3e-12`. Timing is noisy, but the regression persisted. Retain the prior force loop for workloads where Type 2 setup costs dominate.

Other validation completed:

| Command | Result |
|---|---|
| `make -C apps/CCOptimization test-finufft-forces` | Passed, including a fresh rebuild |
| `make -C apps/CCOptimization test-finufft` | All multi-species checks passed, including a fresh rebuild |
| `bash potential/tests/finufft_type1/run_tests.sh` | All three variants passed: 40,977 / 40,999 / 40,999 checks |
| `make -C apps/CCOptimization -j8 all && python3 apps/CCOptimization/tests/test_cli.py` | Build and CLI tests passed |

Inspection supports the force formula, positive transform sign, component layout, duplicate-mode accumulation, cold-cache initialization, exact-type guard, and success-only result writes. No numerical correctness defect was observed. Allocation-failure and nonzero-status fallbacks were inspected but not fault-injected.

**Review limitation:** `git diff c352d19 HEAD` failed because this copy has no Git metadata. I reviewed the listed files; the timing comparator reconstructs the specified prior force loop rather than building `c352d19`. Exact baseline Type 1 arithmetic and default-path performance identity therefore remain unverified.