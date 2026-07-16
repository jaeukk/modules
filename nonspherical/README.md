# nonspherical

Nonspherical-particle shapes and dispersions. Library module (no `main`).

This module **extends core's `SpherePacking` object to nonspherical (polygonal /
polyhedral) packings**: where `SpherePacking` describes a packing of spheres, the
`Dispersion` container here describes a packing of arbitrary `Shape`s (polygons,
polyhedra, …). Correspondingly, it **extends the spectral-density calculation
$\tilde{\chi}_{V}(k)$ from sphere packings to polygonal packings** — the
sphere-packing $\tilde{\chi}_{V}(k)$ lives in `pairstat/SpectralDensity`, and
`SpectralDensity_Dispersion` here is its nonspherical analog, built from each
particle's shape form factor instead of the sphere form factor $\tilde m(k;R)$.

Extracted from the legacy `NonsphericalParticles/` subdir — the core-only pieces
needed by consumers (shapes + dispersions + the dispersion spectral density);
the order-metric / network pieces are left for later.

| File | Purpose |
|------|---------|
| `ParticleShapes.{h,cpp}` | Shape descriptors (spheres, ellipses, polygons, …). |
| `FormFactors.cpp` | `GetFormFactor` implementations for the shapes (PlanarFigure, Polyhedron, …). |
| `Dispersions.{h,cpp}` | `Dispersion` container of nonspherical particles + `ReadDispersion`/I/O. |
| `SpectralDensity_Dispersion.{h,cpp}` | Orientationally-averaged spectral density $\tilde{\chi}_{V}(k)$ of nonspherical dispersions (the nonspherical analog of pairstat's spherical `SpectralDensity`). Depends on `core` only. |

## Dependencies
- **Lower modules:** `core`.
- **External:** GSL (via core).

## Build
```sh
make            # compile objects
make clean
```
