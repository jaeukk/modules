# Round 1 audit: Option B (`finufft-sequential`)

**Verdict: no blocking findings in the C++ evaluation logic.** I found three blocking items outside it: one test that is wrong, one generated-file problem to fix before commit, and missing documentation that the standing rules require. All three are cheap to fix, and none touches the numerics.

**How I checked:** I could only read and search files. I could not run git, compile, or run any test. So I checked that Option A is unchanged by reading the current source, not by diffing against `ab13cfc`. Build warnings and test results rest on Codex's logs.

## Blocking findings

### B1. `tests/threading/test_cli.py:18` – the direct-only branch can never run, so the test fails in the default build
- **Trigger:** build with plain `make` (`USE_FINUFFT=0`), then run `python3 tests/threading/test_cli.py`.
- **Why:** `SetRhoBackend` throws `"...backend not available; FinufftType1 requires compiling with -DCCPOTENTIAL_USE_FINUFFT..."` (`CollectiveCoordinatePotential.h:165`). `main.cpp:1686` catches it, prints `Multicomponent run failed: ...` to stderr and returns 1. That text never contains `"not compiled"`, so the test skips its fallback branch and fails at `assert result.returncode == 0`.
- **Impact:** the test fails in the default build. Its docstring and PASS line also claim "missing-backend behavior" coverage that never runs. Separately, the CLI prints the same message for `finufft` and `finufft-sequential`, so the test cannot tell whether Option B was actually selected.
- **Minimal fix:**
  ```python
  if "backend not available" in result.stderr:
      assert result.returncode != 0
  ```
  Optionally, also print the selected schedule in `main.cpp:1349/1486` (e.g. `"Species schedule: sequential (Option B)"`) and assert it in the test.

### B2. Pre-commit: `tests/threading/__pycache__/` is untracked and not ignored
- **Trigger:** `git add apps/CCOptimization/tests/threading/` would commit `compare.cpython-313.pyc` and `test_cli.cpython-313.pyc`. The root `.gitignore` has no `__pycache__` entry.
- **Fix:** add `__pycache__/` to `.gitignore`, or delete the directory before staging.

### B3. New functions lack the required Doxygen / NumPy docstrings
- `MultiStealthyPotential.h:127–132` (`EvaluateSequentialSpecies`): no `\param` or `\return`.
- `MultiStealthyPotential.h:228–239`: `SetSequentialFinufftSpecies` has no `\param enabled`, and `GetSequentialFinufftSpecies` has no comment at all.
- `MultiStealthyPotential.cpp:806` and `:818`: these definitions have no doc blocks, although every other definition in the file has one.
- The header docs for `Energy` (370–375), `AllForce` (383–389) and `EnergyAndAllForce` (392–399) describe only Option A and are now inaccurate when Option B is on.
- `tests/threading/compare.py`: `choose_cpus` and `main` have no NumPy-style docstrings.
- **Suggested text for `EvaluateSequentialSpecies`:**
  ```cpp
  /** \brief Option B: evaluate species one after another with internally threaded FINUFFT while one
   *  asynchronous worker computes the contact terms (see SetSequentialFinufftSpecies).
   *  \param[out] forces  Per-particle total forces (resized to the particle count), or nullptr for energy only.
   *  \param energy       Compute the total energy; false for AllForce.
   *  \return             The total energy if energy is true, else 0.
   *  \throw std::invalid_argument for a zero thread budget; std::runtime_error for an unassigned particle;
   *         any species or contact exception (the contact worker is joined before it propagates). */
  ```
- For the other three header docs, add one sentence each: "With SetSequentialFinufftSpecies(true), this dispatches to EvaluateSequentialSpecies."

I recommend treating B3 as blocking because it is a standing rule and the fix is comments only, so round 2 stays trivial.

## Checked with no finding
- **Future lifetime and exceptions:** `contactWork` is declared after every buffer it captures, so on an exception it is destroyed (and the `std::async` worker joined) first. `get()` passes on contact errors. `RestoreThreads` resets `ParallelNumThread` even though `OnOneThread` (which predates this change and is not exception-safe) does not, and the ThrowingCC test really exercises this: without the restore the value would be 1, not 4.
- **Data races:** the contact worker only writes `contactEnergy` and `*forces`, which were sized before launch, and nothing else touches them until `get()`. The species loop only reads `pConfig` (coordinates, volume, basis). The cell list's mutable neighbor caches are touched only by the single contact worker, and `PrepareIterateThroughNeighbors(rc)` warms them up beforehand.
- **Deadlocks:** there are no locks or waits, apart from `get()` on a worker that itself waits on nothing.
- **FINUFFT plans and threads:** plans are made serially before the worker launches. Changing `FinufftNumThreads` invalidates the cache fingerprint, which destroys and rebuilds both Type 1 and Type 2 plans. `PrepareFinufftForces(1)` keeps a single batched plan with `ntransf = dim`. Type 1 and Type 2 execution failures and Type 2 plan failures fall back to direct sums.
- **Gradients:** `ContactEnergyAndForce` and `WithSpeciesForce` use the same expressions as before, and the finite-difference test uses the correct relative-coordinate projection.
- **Thread budget and Option A:** five runnable threads on four CPUs is by design and was not flagged. The `OMP_THREAD_LIMIT=1` serial path is covered. Option A's code paths only gained the dispatch lines at the top of each entry point. No caller wraps these calls in an OpenMP region: `main.cpp`, `StructureOptimization.cpp` and `MD_System.cpp` have no active `omp` pragmas.

## Non-blocking notes
1. **Derived species can still run FINUFFT inside an OpenMP region** (`.cpp:881–882`). Take a derived species that inherits base `GetRho` (which any class opting in via `SupportsIndexedView` must), with FINUFFT and more than one Type 1 block. Base `GetRho`'s `omp parallel for if(blocks>1)` then runs FINUFFT with `nthreads=budget` inside a team. Results stay correct; it only breaks the "no enclosing region" contract. Only test-defined derived classes can reach this. Optional fix: set `pot.ParallelNumThread = eligible ? 1 : transformThreads;` before `pot.GetRho()`; `RestoreThreads` undoes it afterwards.
2. **Untested branch:** the multi-block Type 1 loop in Option B is probably never exercised, since 100 particles per species likely yields one block. The risk is low.
3. **Thread creation failure:** if `std::async` cannot start a thread, its `std::system_error` propagates instead of falling back to running contact inline. Acceptable.
4. **Style nit:** `RestoreThreads` has a destructor but implicit copy operations (C++ Core Guidelines C.21). It is local and never copied, so this is harmless.

**Round 2 scope:** B1 (test), B2 (commit hygiene) and B3 (comments/docstrings) only.