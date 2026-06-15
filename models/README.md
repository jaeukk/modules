# models

Generators of point patterns / packings. Library module (no `main`).

The sources here are the **top-level** copies — the legacy `Models/` subdir
`GenerateConfigs` is a 129-line stub, so the 479-line top-level implementation is
canonical (see `../MODULE_SOURCES.md`).

## Contents

| File | Purpose |
|------|---------|
| `GenerateConfigs.{h,cpp}` | Configuration factory: lattices, perturbed lattices, vacancies, RSA, hard-sphere fluids, paracrystals. |
| `RandomSequentialAddition.{h,cpp}` | RSA with voxel-based spatial acceleration. |
| `HardSphereFluids.{h,cpp}` | Monte-Carlo hard-sphere fluids (canonical / isobaric, mono / polydisperse). |
| `RandomOrganization.{h,cpp}` | Berthier active–passive random-organization model. |

## Dependencies

- **Lower modules:** `core`.
- **Optional:** `pairstat` — structure-factor-driven generation, gated behind
  `-DSTRUCTUREFACTOR_INCLUDED`. Enable with `make WITH_SK=1`.
- **External:** GSL, OpenMP.

## Build

```sh
make            # compile objects (core only)
make WITH_SK=1  # also enable StructureFactor-based features (adds pairstat)
make clean
```
