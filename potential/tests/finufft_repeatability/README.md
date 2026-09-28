# FINUFFT repeatability regression

This is a finite regression for `ShiftedCCPotential::SetFinufftNumThreads` /
`GetFinufftNumThreads`. The default is one FINUFFT thread, independent of
`ParallelNumThread`. The tests check two things:

- the API contract;
- that repeating the same optimization with the same backend gives byte-identical
  results on this machine, build, and settings.

They do not claim a general reproducibility guarantee. They also do not test whether
the FINUFFT and Direct backends end at the same particle positions.

## Files

| File | Purpose |
|---|---|
| `run_repeatability.py` | Driver: fresh build, API tests, probe, 12 optimizations, checks, report |
| `finufft_observer.cpp/.h` | Link wrapper for `finufft2d1` (`-Wl,--wrap=finufft2d1`) that only observes. It passes every argument, including the options pointer, through unchanged. It records the call count, each call's `opts->nthreads`, null options, and nonzero statuses, then prints one `FINUFFT_OBSERVER ...` line to stderr at exit |
| `api_test.cpp` | Small API test, built with and without `-DCCPOTENTIAL_USE_FINUFFT` |

The optimizations reuse `../contact_reproducibility/experiment.cpp` unchanged. That file
wraps `../five_component/experiment.cpp` with its defaults: N2000, Species5, Threads4,
and 80 modes (chi_i=80/798). The build uses `../finufft_type1/run_comparison.py::build`,
and the displacement metric comes from `../contact_reproducibility/run.py::compare`.
The driver loads these helpers with `importlib` and does not modify `sys.path`. It does
not reuse any old build strategy or old objects.

## Run

Needs the same toolchain as the earlier tests: `g++`, GSL/NLopt under
`/home/jaeukk/install`, FINUFFT under `<task>/control/native-finufft`, Python 3 with
`numpy`, and at least four physical cores.

```bash
python3 _multiagent/tasks/2026-09-22-graft-finufft/source/modules/potential/tests/finufft_repeatability/run_repeatability.py \
    --output /ABSOLUTE/NEW/DIR/outside/source
```

`--output` must be an absolute path that does not exist yet and lies outside the task's
`source/` tree. Everything runs serially. The exit code is 0 only if every required
check passes, 1 if any fails (including build or harness errors), and 2 for invalid
arguments. The driver never changes existing results or restarts any sweep.

## What it does

1. **Environment.** Uses `five_component/run_study.py::environment`: affinity set to four
   distinct physical cores; `OMP_NUM_THREADS=OMP_THREAD_LIMIT=4`, `OMP_WAIT_POLICY=PASSIVE`,
   `OMP_DYNAMIC=false`, `OMP_PROC_BIND=false`.
2. **Fresh build.** `run_comparison.build` compiles every current source with
   `-DCCPOTENTIAL_USE_FINUFFT`, including `CollectiveCoordinatePotential.h` as it is now
   and the fixed-order contact sum in `MultiStealthyPotential.cpp`. The link is then
   rerun with `finufft_observer.o` and `-Wl,--wrap=finufft2d1`. Nothing wraps or changes
   options: the harness never calls `SetFinufftNumThreads`, so it exercises the
   production default.
3. **API tests** (`api/`). The macro-enabled variant links the fresh core,
   `Potential`, and `MultiStealthyPotential` objects plus the observer and libfinufft
   (no `main` object, no second full rebuild). The macro-disabled variant recompiles
   only `api_test.cpp` and `MultiStealthyPotential.cpp`, without the macro or FINUFFT
   include paths, and links no FINUFFT library. The core objects do not include the
   header, so there is no ODR mismatch. `ldd` must show no libfinufft. The checks cover:
   - default 1 (also for `ShiftedCCPotential_S2` and the species potentials);
   - Direct as the default backend;
   - getter round trip;
   - 0, -1, and `INT_MIN` throw `std::invalid_argument` and leave the option, backend,
     tolerance, max modes, `RhoValid`, and cached rho bitwise unchanged, with no
     FINUFFT recomputation;
   - a valid setter invalidates `RhoValid`, even when the value is unchanged;
   - copy construction, assignment, `clone`, `LogCCPotential`, `ShiftedCCPotential_ForMC`,
     and `_withWeight` copies keep the option, and copies do not share it;
   - `ParallelNumThread` changes, set directly or propagated by `MultiStealthyPotential`,
     leave the option alone;
   - the observer sees `opts.nthreads` 1 by default and 4 when 4 is set explicitly, and
     `ParallelNumThread=4` is not forwarded;
   - with the macro off, FINUFFT cannot be selected and the option is stored but unused.
4. **Probe** (`probe.*`). Runs 100 fixed-coordinate evaluations per backend, with
   independent reference and finite-difference checks inside the harness. The energy
   range and maximum force change must be exactly 0 for both backends. Every observed
   FINUFFT call must request nthreads=1.
5. **Optimizations** (`runs/`). Covers seeds 20260928/29/30 × two separate-process
   repeats × {direct, finufft} = 12 runs, with the backend order alternating as in the
   earlier studies. Initial and final files, modes, direct structure factors, the JSON
   receipt, and stdout/stderr are kept. Each receipt records the selected and actual
   backend, optimizer statuses and evaluation counts, direct-reference acceptance, the
   observer receipt, file hashes, optimization wall time, and process wall time.
6. **Provenance** (`manifest.json`, `build/`, `extra_build/`). Records every compile and
   link command and log, hashes of sources/headers/helpers, binaries, and libfinufft,
   `ldd`, the compiler, platform, affinity, and thread environment.

## Required checks (fail closed)

- Thread environment as above; driver and build completed.
- Both API variants pass; the macro-disabled binary has no libfinufft dependency.
- The probe passes; the 100 repeats are identical for both backends; observed nthreads=1.
- Per run:
  - exit 0;
  - expected N, threads, backend, 80 modes, and chi_i=80/798;
  - direct-reference ground state (E<=1e-11, max S<=1e-10, max force<=1e-7, max
    overlap<=1e-6);
  - no pass cut short by the 600 s budget.
- FINUFFT runs: `finufft_calls == rho_calls > 0` (no fallback); observer calls
  >= harness calls; every call requested nthreads=1 and returned 0. Direct runs: zero
  FINUFFT calls.
- Per seed:
  - all four runs have byte-identical `initial.tsv` and `modes.tsv`;
  - within each backend, the two repeats have byte-identical `final.tsv`, identical pass
    records (status, evaluations, energy), identical rho/FINUFFT counters, and identical
    observer receipts.
- Tracked sources are unchanged between the start and end of the run.

Descriptive only: cross-backend displacement (fixed labels, periodic, translation-aligned),
same-backend displacement, optimizer success codes (1–4), and timings.

Results: `summary.json` (all checks with details), `REPORT.md`, and the directories above.

## Assumptions and limits

- The FINUFFT API is taken from the agreed interface, which the header already contained
  when these tests were written: invalid values throw `std::invalid_argument`, and every
  valid call invalidates the cache.
- Direct-reference acceptance and the absence of time-budget truncation are required.
  A timed-out run depends on wall-clock time, so comparing its repeats would be
  meaningless. Optimizer codes are recorded, not required.
- Only `finufft2d1` is observed, since the harness is two-dimensional.
- A passing run supports repeatability for these 12 cases on this machine and build.
  It says nothing about other thread counts, sizes, platforms, or the FINUFFT/Direct
  agreement of endpoints.
