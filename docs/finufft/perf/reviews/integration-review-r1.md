**Two BLOCKING findings and one FOLLOW-UP.**

1. **BLOCKING — Copying an evaluated potential can produce incorrect Type 2 forces.**

   **Location:** `potential/CollectiveCoordinatePotential.h:431`, `:437`, `:1008`.

   **Trigger:** Evaluate a FINUFFT `ShiftedCCPotential`, copy it, then call `AllForce()` on the copy without rebinding its unchanged configuration.

   The implicit potential copy preserves `RhoValid=true` and the FINUFFT backend marker, while `FinufftPlanCache` starts empty. Preparation allocates zero-filled angles, but `TryFinufftAllForce()` skips `GetRho()` because rho remains valid. Type 2 therefore evaluates the copied rho at incorrect particle positions. Assignment has the same validity mismatch.

   **Evidence:** I compiled a reproduction using the existing `test_finufft_forces.cpp` helpers, with a skew 2D cell, N=2000, χ=0.4, and four threads. The essential sequence was:

   ```cpp
   a.AllForce(fa);
   ShiftedCCPotential b(a);
   b.AllForce(fb);
   ```

   Running `/tmp/multiagent-review-3qarxbdw/tmp/review_repro` produced:

   ```text
   COPY gate=1 type2=1 max_error=399.879 scale=404.873
   COPY_VS_DIRECT max_error=399.879 scale=404.873
   ```

   This exceeds the force-error threshold substantially. The existing copy test at `apps/CCOptimization/tests/test_finufft_plan_cache.cpp:525` manually clears `RhoValid`, masking the defect.

   **Suggested correction:** Invalidate copied rho whenever its associated FINUFFT angle cache is discarded.

2. **BLOCKING — The contact-energy schedule slows the default Direct backend on small systems.**

   **Location:** `apps/CCOptimization/MultiStealthyPotential.cpp:392`.

   **Trigger:** The existing benchmark’s ordinary random packing with d=3, N=100, S=5, χₛ=0.4, and T=4 or 8.

   Replacing static scheduling with `schedule(dynamic, 64)` gives 100 particles only two chunks. Consequently, at most two workers compute contact energy, including when every species uses Direct.

   **Evidence:** I reconstructed the pre-integration implementation by reversing the supplied diff for the potential header and wrapper files, then built the unchanged benchmark against it. Repeated commands were:

   ```bash
   OMP_PROC_BIND=close OMP_PLACES=cores /tmp/multiagent-review-3qarxbdw/tmp/bench_base eval 3 100 5 .4 4 1
   OMP_PROC_BIND=close OMP_PLACES=cores apps/CCOptimization/bench_finufft.out eval 3 100 5 .4 4 1
   ```

   I repeated these with T=8. Direct-backend median evaluation times, milliseconds:

   | Threads | Pre-integration, three runs | Integrated, three runs |
   |---|---|---|
   | 4 | 0.1315, 0.1308, 0.1231 | 0.1716, 0.1522, 0.1815 |
   | 8 | 0.1235, 0.1234, 0.1300 | 0.1655, 0.1740, 0.1652 |

   A diagnostic build changing only that schedule back to static removed the regression. In another interleaved T=8 comparison, pre-integration measured 0.1168–0.1415 ms, integrated 0.1567–0.1703 ms, and static 0.0848–0.1019 ms.

   **Suggested correction:** Preserve static contact scheduling for the Direct path.

   **Baseline limitation:** Git metadata was unavailable; this comparison uses the supplied diff’s `dd2fb45` source state, not an independently verified checkout of `c352d19`.

3. **FOLLOW-UP — Destroying a successful Type 2 plan prevents its later recreation.**

   **Location:** `potential/CollectiveCoordinatePotential.h:378`.

   **Trigger:** Evaluate with FINUFFT, switch to Direct and evaluate, then switch back to FINUFFT with the same configuration size, constraints, and settings.

   Destruction clears `plan2` but leaves `status2=0`. The recreation condition requires `FinufftNotCalled`, so forces remain Direct until a fingerprint change rebuilds the cache.

   The reproduction reported:

   ```text
   BACKEND_TOGGLE gate=1 type2=0
   REACCEPT gate=1 prepared=0
   ```

   This loses acceleration and makes backend selection depend on evaluation history. I did not establish a slowdown against the original FINUFFT baseline specifically for this transition. Reset the successful status when deliberately destroying the plan, while retaining actual failure caching.

**Validation performed**

- `make -C apps/CCOptimization test-finufft test-finufft-forces test-species-concurrent test-indexed-view test-finufft-plan`: multi-species oracle, 163 cache checks, and all three force-test variants passed. The command stopped at the concurrency test’s thread-count assertion: `peak of 10 OS threads exceeds the 8-thread team`; its other 2,257 checks passed.
- Separate indexed-view run passed. Default concurrency rerun, `OMP_WAIT_POLICY=passive`, and `OMP_MAX_ACTIVE_LEVELS=4` each passed all 2,258 checks.
- `bash potential/tests/finufft_type1/run_tests.sh`: all three variants passed.
- Header-prescribed ASan build: LeakSanitizer could not run under ptrace. With `ASAN_OPTIONS=detect_leaks=0`, all 163 checks passed without reported address errors. Leak freedom remains unverified.
- `make -C apps/CCOptimization -j8 all && python3 apps/CCOptimization/tests/test_cli.py`: passed.
- FINUFFT benchmark build and `bench_finufft.py --quick`: passed.
- Extended benchmark matrix covered d=2/3, N=100/2000, S=5, χₛ=0.1/0.4, T=4/8. Integrated FINUFFT improved over reconstructed pre-integration FINUFFT throughout; the reproducible Direct regression is identified above.