# nonspherical

Nonspherical-particle shapes and dispersions. Library module (no `main`).
Extracted from the legacy `NonsphericalParticles/` subdir — only the core-only
pieces needed by consumers (shapes + dispersions); the spectral-density / form-
factor / order-metric / network pieces are left for later.

| File | Purpose |
|------|---------|
| `ParticleShapes.{h,cpp}` | Shape descriptors (spheres, ellipses, polygons, …). |
| `FormFactors.cpp` | `GetFormFactor` implementations for the shapes (PlanarFigure, Polyhedron, …). |
| `Dispersions.{h,cpp}` | `Dispersion` container of nonspherical particles + `ReadDispersion`/I/O. |

## Dependencies
- **Lower modules:** `core`.
- **External:** GSL (via core).

## Build
```sh
make            # compile objects
make clean
```
