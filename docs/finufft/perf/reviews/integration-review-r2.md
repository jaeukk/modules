**No blocking findings under the agreed baseline rule. One FOLLOW-UP regression.** Review scope was limited to `integration-r2-fix.diff`.

1. **FOLLOW-UP — Larger concurrent contact chunks slow large-N Energy relative to round 1.**

   **Location:** `apps/CCOptimization/MultiStealthyPotential.cpp:409`, `:417`.

   **Trigger:** Existing benchmark setup: random 3D packing, N=10,000, S=5, χₛ=0.4, FINUFFT, four or eight threads.

   **Impact:** Repeated Energy measurements were slower than round 1. Changing only the concurrent schedule back to `dynamic,64` reduced the slowdown.

   I reconstructed round 1 by reversing the fix diff, built identical benchmark harnesses, and ran:
   ```bash
   python3 review/fasttiming.py
   ```
   This alternated execution order across five runs, with 120 measured evaluations per process and `OMP_PROC_BIND=close OMP_PLACES=cores`.

   Median-of-run-medians, milliseconds:

   | N | Threads | Round 1 | Fixed | Fixed with chunk 64 |
   |---|---:|---:|---:|---:|
   | 100 | 4 | 0.226 | 0.228 | 0.238 |
   | 100 | 8 | 0.244 | 0.246 | 0.244 |
   | 10,000 | 4 | 6.927 | 7.328 | 7.090 |
   | 10,000 | 8 | 6.686 | 7.869 | 6.392 |

   The large cases regress approximately **6% and 18%** against round 1. Small-case differences are within timing noise. Prefer retaining the previous concurrent chunk size while keeping the restored static Direct schedule.

   This is FOLLOW-UP because FINUFFT remains faster than the reconstructed pre-integration baseline; it does not meet the agreed blocking criterion requiring a slowdown against that baseline.

**Fix verification**

- Evaluated-object copy, clone, and assignment preserved source force bits at four and eight threads in an added skew-cell N=2,000 reproduction. Direct comparisons passed.
- Injected Type 1 fallback, subsequent recovery, reversed indexed-view rebinding, invalidated constraint weight/mode edits, and tolerance changes passed. Maximum observed force error was approximately `7.4e-11` on force scale `74.2`.
- Backend toggling restored Type 2. Injected Type 2 makeplan failure remained cached across decline/accept and backend toggles.
- Direct Energy matched reconstructed baseline bits at N=100 and N=10,000 with four/eight threads. The round-1 small-N slowdown disappeared. Timings were broadly comparable; shared-machine noise prevents claiming exact performance identity.
- The revised sampler rejected a deliberate 12-thread team against a four-thread limit:
  ```text
  FAIL: deliberately oversized team: peak=12
  failures=1
  ```

**Commands and results**

- `make -C apps/CCOptimization test-finufft test-finufft-forces test-species-concurrent test-indexed-view`: passed, including 168 plan-cache checks and all force-test variants.
- Three additional concurrency runs: passed. Across all six runs, each reported 2,337 checks passed and peak eight threads.
- `bash potential/tests/finufft_type1/run_tests.sh`: all three variants passed.
- Header-prescribed ASan build: LeakSanitizer failed because ptrace is unsupported. Rerunning with `ASAN_OPTIONS=detect_leaks=0` passed all 168 checks without address errors; leak freedom remains unverified.
- `make -C apps/CCOptimization -j8 all && python3 apps/CCOptimization/tests/test_cli.py`: passed.
- Added `review/extra` checks passed; `review/peak` failed as intended. `review/timing.py` and `review/fasttiming.py` completed the baseline and isolated-schedule comparisons.

Git metadata was unavailable. Baseline comparisons use sources reconstructed by reversing the supplied diffs, not independently verified commit checkouts.