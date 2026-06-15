# pairstat

Pair / two-point statistics of point patterns and packings. This is the **lean,
decoupled** version of these routines (from the legacy `PairStat/` subdir), with
no dependency on `NonsphericalParticles`.

> If you need the nonspherical / spectral-density extensions that the legacy
> top-level `StructureFactor.cpp` / `NumberVariance.cpp` carried, see
> `../MODULE_SOURCES.md` — those features were intentionally dropped here and
> should be ported behind a `USE_SHAPES` guard if a consumer requires them.

## Contents

| File | Purpose |
|------|---------|
| `PairCorrelation.{h,cpp}` | Radial distribution function g(r), pair-correlation histograms. |
| `StructureFactor.{h,cpp}` | Static structure factor S(k), isotropic & directional. |
| `NumberVariance.{h,cpp}` | Local number variance σ²(R) (uses Boost gamma). |
| `Computation.{h,cpp}` | High-level driver tying the above together for the CLI. |
| `ContactStatistics.{h,cpp}` | Contact-number / nearest-neighbour statistics. |

## Dependencies

- **Lower modules:** `core`.
- **External:** GSL, Boost (header-only special functions), OpenMP.

## Build

```sh
make            # builds module objects + pairstat.out (CLI)
make objs       # object files only
./pairstat.out  # dispatches on a command argument (prints usage otherwise)
```

The makefile builds the `core` module's objects automatically as a prerequisite.
