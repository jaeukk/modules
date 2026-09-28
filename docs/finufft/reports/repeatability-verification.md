## PASS — finufft-repeatability-fix

Own executions (stdout-only, affinity `0,2,4,6`, requested OpenMP environment):

```bash
env OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 OMP_DYNAMIC=false OMP_WAIT_POLICY=PASSIVE OMP_PROC_BIND=false taskset -c 0,2,4,6 .../api_test_finufft
env OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 OMP_DYNAMIC=false OMP_WAIT_POLICY=PASSIVE OMP_PROC_BIND=false taskset -c 0,2,4,6 .../api_test_direct
env OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 OMP_DYNAMIC=false OMP_WAIT_POLICY=PASSIVE OMP_PROC_BIND=false taskset -c 0,2,4,6 .../experiment_observed probe
```

Results: enabled API passed 75 checks; disabled API passed 50. Probe passed: both backends had exactly zero energy/force range over 100 repeats. Its observer recorded 565 FINUFFT calls, all `nthreads=1`, zero null options/status failures. The enabled API’s `1:5,2:2,3:1,4:1` receipt is expected from explicit setter tests, confirming values propagate; it is not the default probe path. Disabled binary has no `libfinufft` dependency.

Source review: `FinufftNumThreads` initializes to 1, validates `nthreads >= 1`, invalidates only rho cache, and is passed directly to `opts.nthreads`. It remains distinct from `ParallelNumThread`; probe/runs retain four threads elsewhere. The observer forwards the original options pointer unchanged. The fixed-order contact sum remains per-particle slots followed by ordered serial accumulation; no diff versus the task source baseline.

Inspected conductor evidence:

- All 12 receipts present; all current source hashes match manifest hashes.
- For every seed, all four initial/mode files are byte-identical; same-backend finals, pass records, rho/FINUFFT counters, and observer receipts are byte-identical.
- FINUFFT runs have `finufft_calls == rho_calls > 0`, all observer calls successful and `nthreads=1`; Direct runs have zero FINUFFT calls.
- Independently recalculated acceptance criteria from receipts: all runs satisfy energy ≤1e-11, max S ≤1e-10, max force ≤1e-7, overlap ≤1e-6; no status 6 or 600-second truncation.
- Optimization and process timings are separately recorded for every run.

No cross-backend endpoint equality was required or assumed.

Service status: `systemctl` is inaccessible in this sandbox (“Failed to connect to bus”); process inspection showed no matching active sweep/FINUFFT-repeatability process beyond the inspection command itself.