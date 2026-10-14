**No blocking findings.** No new FOLLOW-UP findings within the fix’s scope.

The supplied diff preserves the mapping checks, arithmetic order, grid limits, and transform indexing. `FinufftEligible()` shares those checks with the transform path. In the exercised cases, its prediction matched the backend actually used.

Fallback species now retain their full OpenMP team in both `Energy()` and cold `AllForce()`. All-fallback results were bitwise identical to Direct. Mixed eligible/fallback cases remained repeatable at T=4 and T=8.

Validation completed:

| Command | Result |
|---|---|
| `make -C apps/CCOptimization test-species-concurrent` | 2,258 checks passed under each of default waiting, passive waiting, and nesting enabled. |
| `make -C apps/CCOptimization test-finufft` | All multi-species checks passed. |
| `bash potential/tests/finufft_type1/run_tests.sh` | All three variants passed: 40,960 / 40,982 / 40,982 checks. |
| `make -C apps/CCOptimization -f review_r2.mk USE_FINUFFT=1 review_hash.out` followed by two runs of `OMP_WAIT_POLICY=passive apps/CCOptimization/review_hash.out` | Both runs passed all 2,258 checks. All four mixed-case energy/force hashes matched across processes; they also matched between T=4 and T=8. |

I reproduced the round-1 performance case with N=4,000, S=1, T=8, χ=0.4, shear=100, and contact disabled. The disposable harness invalidated rho before each timed evaluation and compared the fixed wrapper, Direct wrapper, and a serial species call using FINUFFT selection with the full team:

```text
make -C apps/CCOptimization -f review_r2.mk USE_FINUFFT=1 review_r2.out
OMP_WAIT_POLICY=passive apps/CCOptimization/review_r2.out

fallback_direct_bitwise=1 eligible=0 modes=3196
```

Four batches of 12 evaluations produced these median times:

| Batch | Fixed fallback | Direct | Serial species control |
|---|---:|---:|---:|
| 1 | 71.40 ms | 69.65 ms | 67.88 ms |
| 2 | 74.29 ms | 69.69 ms | 64.74 ms |
| 3 | 60.38 ms | 61.93 ms | 61.67 ms |
| 4 | 61.13 ms | 58.83 ms | 61.43 ms |

The previous approximately 5× fallback slowdown was not reproduced; performance is now comparable to the full-team controls. Small timing differences remain inconclusive on this shared machine.

Git metadata is absent: `git status` and `git log` returned “not a git repository.” Consequently, historical bitwise identity with c352d19 could not be directly verified. The serial species control reconstructs its dispatch behavior for this performance case; it is not a historical checkout. Runtime FINUFFT-error fallback was not fault-injected and remains the documented exception.