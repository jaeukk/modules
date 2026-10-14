**Request changes: one BLOCKING performance regression.** Numerical correctness and repeatability checks passed.

Git metadata was unavailable (`git diff c352d19 HEAD` returned `Could not access 'c352d19'`). I reviewed the listed files and used a serial-dispatch control for comparison; I could not verify historical Direct bitwise identity.

1. **BLOCKING — FINUFFT fallback loses parallel rho evaluation**

   **Location:** `apps/CCOptimization/MultiStealthyPotential.cpp:330` and `:342`.

   **Trigger:** FINUFFT selected, one species, N=4,000, T=8, χ=0.4, contact disabled, and a unit-density periodic cell with basis `(L,0), (100L,L)`, L=√4000. This is a valid skew representation of a square lattice. Its 3,196 constraints trigger the existing sparse-grid fallback to Direct.

   **Impact:** `SpeciesConcurrent()` checks the requested backend; `SpeciesEnergyOnOneThread()` then forces the Direct fallback onto one thread. Energy evaluation slowed approximately **5×**.

   **Evidence:** Built a disposable harness using the new test’s configuration/potential helpers, invalidating rho before every timed evaluation. Ran:
   ```text
   OMP_WAIT_POLICY=passive ./review_fallback.out

   N=4000 S=1 T=8 chi=.4 shear=100 modes=3196 actual_backend=Direct
   serial-species baseline rho/energy median_ms=50.1427
   option-C median_ms=264.672003
   ```
   Four batches of 12 evaluations gave **49.48–52.97 ms** for the prior per-species call versus **262.91–266.24 ms** for option C.

   I also compiled the same wrapper with only `SpeciesConcurrent()` changed to return false:
   ```text
   OMP_WAIT_POLICY=passive ./review_serial.out
   ```
   Full wrapper energy medians were **48.38–51.81 ms**, confirming the dispatch change causes the regression. This control reconstructs the serial dispatch; it is not a checkout of c352d19.

   The acknowledged fallback limitation meets the agreed slowdown blocking bar. Preserve parallel Direct evaluation when fallback occurs.

2. **FOLLOW-UP — Thread-count test masks normal DUCC0 pool creation**

   **Location:** `apps/CCOptimization/tests/test_species_concurrent.cpp:246`.

   **Trigger:** Run without `DUCC0_NUM_THREADS` or `OMP_NUM_THREADS` set. The test itself forces `DUCC0_NUM_THREADS=1`; the application does not.

   **Impact:** The claimed peak OS thread count ≤ T is only demonstrated with that extra environment setting.

   **Evidence:** Built a test copy removing that assignment and adding numerical-result hashes:
   ```text
   env -u DUCC0_NUM_THREADS -u OMP_NUM_THREADS \
     OMP_WAIT_POLICY=passive ./review_uncapped.out

   FAIL: peak of 15 OS threads exceeds the 8-thread team
   peak OS threads (excluding sampler): 15
   FAILED: 1 of 1378 species-concurrent checks failed
   ```
   All numerical checks passed. This is **FOLLOW-UP**, because the extra idle pool alone does not establish a qualifying crash, numerical defect, or slowdown.

Other validation performed:

- `make -C apps/CCOptimization test-species-concurrent`: **1,378 checks passed** in each of default, passive-wait, and nesting-enabled modes.
- `make -C apps/CCOptimization test-finufft`: passed.
- `bash potential/tests/finufft_type1/run_tests.sh`: all three variants passed—40,960 / 40,982 / 40,982 checks.
- `make -C apps/CCOptimization -j8 all` and `python3 apps/CCOptimization/tests/test_cli.py`: passed.
- Two fresh processes running the hash-instrumented test with `DUCC0_NUM_THREADS=1`: all 24 energy/force hashes matched; `diff -u` produced no differences.
- `ldd` and symbol inspection of installed `libfinufft.so`: DUCC0 symbols present, no dynamically linked FFTW dependency. No concurrent-call failure observed.

The tests support cold `AllForce` correctness, fresh-potential repeatability, and FINUFFT bitwise equality at T=4 versus T=8 for the exercised cases.