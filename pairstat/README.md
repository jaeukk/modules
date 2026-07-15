# pairstat

Pair / two-point statistics of point patterns and packings. This is the **lean,
decoupled** version of these routines (from the legacy `PairStat/` subdir), with
no dependency on `NonsphericalParticles`. It is a **library-only** module — the
high-level `Computation` driver and its CLI now live in
[`../apps/PairStatCLI/`](../apps/PairStatCLI/README.md).

> The *spherical* spectral-density and volume-fraction-variance routines that the
> legacy top-level `StructureFactor.cpp` / `NumberVariance.cpp` carried have been
> ported back (now in their own `SpectralDensity.{h,cpp}` /
> `VolumeFractionVariance.{h,cpp}`). The *nonspherical* / `USE_SHAPES` extensions
> are still intentionally dropped here — see `../MODULE_SOURCES.md`.

## Contents

| File | Purpose |
|------|---------|
| `PairCorrelation.{h,cpp}` | Radial distribution function g(r), pair-correlation histograms. |
| `StructureFactor.{h,cpp}` | Static structure factor S(k) (isotropic, weighted, directional) and the sphere form factor `Mtilde`. |
| `SpectralDensity.{h,cpp}` | Spectral density χ_V(k) of sphere packings (`SpectralDensity[_Lstorage]`, `IsotropicSpectralDensity[_Lstorage]`). Split out of `StructureFactor`; built from `Mtilde`. |
| `NumberVariance.{h,cpp}` | Local number variance σ²(R) (uses Boost gamma). |
| `VolumeFractionVariance.{h,cpp}` | Local volume-fraction variance σ_V²(R) of sphere packings (the `MC_vfv` family). Split out of `NumberVariance`. |
| `ContactStatistics.{h,cpp}` | Contact-number / nearest-neighbour statistics. |
| `LinealPathFunction.{h,cpp}` | Lineal-path function L(z) and chord-length distribution of sphere packings / two-phase media. |

### Why the splits

`SpectralDensity` is the packing-volume analog of `StructureFactor` (it reuses
the form factor `Mtilde`, which stays in `StructureFactor.h`), and
`VolumeFractionVariance` is the volume analog of `NumberVariance`. Keeping each
quantity in its own translation unit lets a consumer that only needs, say, S(k)
avoid pulling in the sphere-packing spectral-density / volume code. The small
`SkBin` binning accumulator is kept file-local in both
`StructureFactor.cpp` and `SpectralDensity.cpp`.

## Dependencies

- **Lower modules:** `core`.
- **External:** GSL, Boost (header-only special functions), OpenMP.

## Build

```sh
make            # or: make objs — compiles the module object files
make clean
```

This module no longer links an executable; the makefile builds only the object
files, against which consumers (e.g. `../apps/PairStatCLI/`) link. The `core`
module's objects are built automatically as a prerequisite.
