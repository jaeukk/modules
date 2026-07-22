# models

Generators of point patterns / packings. Library module (no `main`).

The code is the reviewed (canonical) module implementation, organized into
per-model subdirectories. The compiled objects keep flat basenames in `obj/`
(`GenerateConfigs.o`, `RandomSequentialAddition.o`, `HardSphereFluids.o`,
`RandomOrganization.o`) and headers resolve flat (`#include <RSA_Gen.h>`,
`<HardSphereFluids.h>`, `<Lattice.h>`, `<GenPerturb.h>`, `<HSF_Gen.h>`,
`<RandomOrganization.h>`, `<GenerateConfigs.h>`) via `$(models_INC)` in
`config.mk`, so no consumer needs subdir-aware includes or makefile edits.

## Contents

| Path | Purpose |
|------|---------|
| `GenerateConfigs.{h,cpp}` | ConfigGen base class, `PoissonPtGen`, `Jammed_Gen` stub, predefined-lattice helpers (`GetPredefinedBasis` / `GetPredefinedCrystals`), `Echo`, `numK_points`. Re-includes the split headers below as a backward-compatible facade. |
| `Lattices/Lattice.{h,cpp}` | `LatticeGen` — periodic lattice packings (impl `#include`d into `GenerateConfigs.o` so consumers linking only that object keep resolving `LatticeGen`). |
| `Pert/GenPerturb.h` | `GenPerturb` (stochastic displacements) + `GenPtsDefect` (vacancies / interstitials); header-only. |
| `RSA/RandomSequentialAddition.{h,cpp}` | RSA with voxel-based spatial acceleration. |
| `RSA/RSA_Gen.h` | `RSA_Gen` ConfigGen wrapper + `PhiC(d)` saturation table. |
| `HSF/HardSphereFluids.{h,cpp}` | Monte-Carlo hard-sphere fluids (canonical / isobaric, mono / polydisperse). |
| `HSF/HSF_Gen.h` | `HSF_Gen` ConfigGen wrapper around `HSF`. |
| `RandOrg/RandomOrganization.{h,cpp}` | Berthier active–passive random-organization model. |

A migrated CLI consumer that drives these generators lives in `../apps/Models`.

## Model notes & references

**Random sequential addition (RSA)** — hard hyperspheres are inserted one at a
time at uniformly random positions, and an insertion is accepted only if it does
not overlap any previously placed sphere; the process continues until no further
sphere can be added (saturation). `RandomSequentialAddition.{h,cpp}` is an
original implementation of the voxel-splitting algorithm of

> 1. Zhang, G. & Torquato, S. Precise algorithm to generate random sequential
>    addition of hard hyperspheres at saturation. *Phys. Rev. E* **88**, 053312
>    (2013),

which guarantees that the RSA packing reaches saturation in a finite time.

**Random organization** — `RandomOrganization.{h,cpp}` implements the
periodically driven random-organization model published in

> 1. Tjhung, E. & Berthier, L. Hyperuniform density fluctuations and diverging
>    dynamic correlations in periodically driven colloidal suspensions.
>    *Phys. Rev. Lett.* **114**, 148301 (2015).

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
