# FINUFFT Type 1 backend for `ShiftedCCPotential::GetRho`

`ShiftedCCPotential::GetRho` in `potential/CollectiveCoordinatePotential.h` can optionally evaluate

rho(k) = sum_j exp(+i k . r_j)

with one FINUFFT Type 1 transform (positive sign, unit strengths, no normalization) instead of the
O(N * Nk) direct sum. The constraint weights `constraints[i].V` are not transform strengths.
`Energy`, `Force`, and `AllForce` are unchanged and consume the cached rho.

## API

```cpp
typedef ShiftedCCPotential::RhoBackend Backend;   // Direct (default) or FinufftType1
ShiftedCCPotential::IsRhoBackendAvailable(Backend::FinufftType1); // true only if compiled in
pot.SetRhoBackend(Backend::FinufftType1);  // throws std::invalid_argument if not compiled in
pot.SetFinufftTolerance(1e-12);            // default 1e-12; must be > 0 and finite
pot.SetFinufftMaxModes(std::int64_t(1)<<22); // memory cap on prod_d (2K_d+1); default 2^22
pot.ParallelNumThread = 2;                 // also the FINUFFT thread cap (opts.nthreads)
pot.GetLastRhoBackendUsed();               // what actually produced the cached rho
pot.GetLastFinufftStatus();                // 0, a nonzero FINUFFT code, or FinufftNotCalled (-1)
```

Every setter invalidates the rho cache (`RhoValid=false`); a setter that throws changes nothing.

## Build

Default builds need no FINUFFT header or library. To enable the backend, define the macro in
**every** translation unit that includes the header (mixing it violates the ODR):

FINUFFT 2.5.1 is installed as a shared library under `$HOME/install`, next to GSL. Its
`finufft_common/defines.h`, which upstream's install omits, was copied there from the source tree.

```sh
g++ -std=c++17 -O2 -fopenmp -DCCPOTENTIAL_USE_FINUFFT \
    -Icore -Ipotential -I"$HOME/install/include" my_app.cpp <core/potential objects> \
    -L"$HOME/install/lib" -Wl,-rpath,"$HOME/install/lib" -lfinufft -lgsl -lgslcblas
```

## Method and fallbacks

For each constraint the header computes n_d = k . a_d / (2 pi) against the current basis vectors
a_d (skew cells included), so k.r_j = sum_d n_d (2 pi s_jd) with fractional coordinates s_j. The
mapping is recomputed on every evaluation because `constraints` is public and mutable. Modes may be
unordered, of either sign, zero, or duplicated. `HAVETWOFIXEDPARTICLES` still adds +2 to Re rho.

`GetRho` falls back to the unchanged direct sum (bitwise-identical result, reported as `Direct`) when:

- the dimension is not 1, 2, or 3;
- there are no constraints or no particles;
- snapping some n_d to an integer could exceed the phase budget: both |n_d - round(n_d)| and the
  roundoff estimate 4*DBL_EPSILON*(1 + sum_c |k_c a_dc| / (2 pi)) must be at most
  min(1e-8, eps) / (8 pi dim), so the snapping phase error stays below min(1e-8, eps)/4;
- |n_d| > 2^30, the grid exceeds `SetFinufftMaxModes`, or the grid exceeds
  `FinufftSparseGridFactor` (64) times (Nk + N) (pathologically sparse);
- an allocation fails, or FINUFFT returns any nonzero status (including warning 1).

The phase budget is deliberately conservative. At the default eps=1e-12 it accepts only modes
with sum_c |k_c a_dc| / (2 pi) up to about 43 (1D), 21 (2D), or 14 (3D); orthogonal cells reach
|n_d| up to those values, skew cells less. For eps below about 2.2e-14*dim every mode is rejected
before FINUFFT is called. Loosening eps (e.g. 1e-10) raises the limit by the same factor.
All mapping checks run before the dense mode array is allocated.

## Scope

Only `ShiftedCCPotential::GetRho` honors the setting, so `LogCCPotential`, `ShiftedCCPotential_S2`,
the `_varBox` variants, and `ShiftedCCPotential_ForMC` inherit it. `ShiftedCCPotential_withWeight` and
`ShiftedCCPotential_S2_MultiConfig` override `GetRho` and keep their direct sums; they copy but
ignore the setting. `ShiftedCCPotential_ForMC::TryMove` still updates rho incrementally with the
direct formula.

## Tests

```sh
bash potential/tests/finufft_type1/run_tests.sh      # run from modules/, or use any path
```

The runner compiles `test_finufft_type1.cpp` in a `mktemp -d` directory as three variants
(direct-only without any FINUFFT flags, FINUFFT, FINUFFT + `HAVETWOFIXEDPARTICLES`), links the
prebuilt `core/obj` and `potential/obj/Potential.o` (build them with `make -C core objs` and
`make -C potential objs` if missing; they must be built without `-DNDEBUG`), runs each with
`OMP_NUM_THREADS=2`, and removes only the files it created. Set `KEEP_WORK=1` to keep them, and
`FINUFFT_PREFIX`, `GSL_PREFIX`, or `CXX` to override the defaults. Because FINUFFT and GSL share
`$HOME/install/include`, the runner also checks that the direct-only build never includes `finufft.h`.

Coverage: default/disabled backend selection and setter validation; analytic sign checks; empty
constraint and particle lists; rho, Energy, Force, and AllForce Direct vs FINUFFT in 1D-3D with
orthogonal, skew, and rotated cells, shuffled modes, and duplicates; fallbacks (non-lattice k,
near-lattice k, a 1D n=1e6+1e-8 mode with N=32000 that only the phase budget rejects, sparse grid,
index overflow, mode cap, 4D, eps=1e-20), each required to skip FINUFFT (`FinufftNotCalled`) and
match Direct bitwise; cache invalidation and mutated constraints; a central finite-difference force
check; derived classes; and a fixed-seed N=500 2D smoke comparison with the 299 shortest +/-k pairs,
chi = 299/(2*(500-1)) ~ 0.300 (timings are informational only). The FINUFFT nonzero-status fallback
path is no longer reached by any test. Tolerances are atol + rtol * scale, with scale
N for rho, sum|V| N^2 / vol for energy, and 4 sum|V| N kmax / vol for forces. Each FINUFFT case
asserts that FINUFFT really ran (`GetLastRhoBackendUsed`, status 0).
