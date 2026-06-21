# Modularization — rules & progress log

This file is the **single source of truth** for the rules agreed during the
modularization of the legacy C++ physics library, and a running log of migration
progress. Keep it updated as work continues.

---

## Rules & conventions

### R1. Non-destructive
The legacy tree (`../c++/` — top-level files and the original consumer subdirs)
is **never modified**. All new work lives under `30_Codes/modules/`. Legacy stays
as a fallback until every consumer is migrated.

### R2. Canonical-copy policy
Prefer the leaner *subdir* copies as canonical (`cores/` → core,
`PairStat/` → pairstat, `Potential/` → potential). Exceptions, where the subdir
copy is a stub or lacks needed features, are recorded in `MODULE_SOURCES.md`.
Notably:
- `GenerateConfigs` is taken from **top-level** (the `Models/` copy is a stub).
- The lean `pairstat` dropped nonspherical / spectral-density / `GetKs` /
  `Mtilde`. A consumer that needs those uses the **full** legacy `StructureFactor`
  (documented per-consumer) or the features get ported behind a flag.

### R3. Build wiring
- Every makefile does `include <relpath>/config.mk` and uses the module path
  vars (`$(core)`, `$(pairstat)`, …). **Never** hardcode `/home/jaeukk/...` or use
  `../File.cpp` to reach the library.
- `config.mk` auto-detects the GSL/NLopt prefix from `$(HOME)` and derives all
  paths from its own location, so the whole `modules/` tree is relocatable.
- Optional heavy features are command-line toggles: `USE_OPENMP` (default 1),
  `USE_CGAL`, `USE_VORO`, `USE_TJ`, plus per-module flags like `WITH_SK`.

### R4. Include style in sources
- Library headers are included with angle brackets: `#include <GeometryVector.h>`
  (resolved by `-I$(module)` paths), **not** `#include "../GeometryVector.h"`.
- A file's own / sibling headers keep quotes: `#include "FFH.h"`.

### R5. Library globals
The `core` module owns the global `Verbosity` (and the other `etc.cpp` globals).
Consumers and their `main` must **not** redefine them.

### R6. Layering
A module may only `-I` modules in a lower or equal tier (see `README.md`).
Library modules live in `modules/<name>/`; migrated application drivers live in
`modules/apps/<name>/`.

### R7. Migration procedure (per consumer)
1. Copy the consumer's own sources into `modules/apps/<name>/`.
2. Rewrite `#include "../X.h"` → `#include <X.h>` (keep sibling quotes).
3. Write a makefile: `include ../../config.mk`, set `-I` for the modules used,
   link the module objects (`make -C $(core) objs`, etc.).
4. Build and run. If it can't build (missing lib / unported feature / unbuilt
   module), migrate the wiring anyway and record it as **blocked** below.

---

## Synced to ../new_c++ (authoritative newer source)

All source files that exist in both `modules/` and the parallel `new_c++/` tree
were compared; the **25 differing files** were replaced with the `new_c++`
versions (newer core/pairstat/potential + app mains). Reconciliations applied so
the module tree still builds:

- **Include wiring re-normalized**: some `new_c++` files use `../X.h`; converted
  back to `<X.h>` to match the flat-module `-I` layout (R4).
- **Verbosity convention adopted from new_c++**: `core/etc.cpp` does *not* define
  it; each executable's `main` does. Re-added `size_t Verbosity` to my own
  drivers (`cli.cpp`, `g2_sk_demo.cpp`, `CHI/main.cpp`, `FirstPassageTime/main.cpp`).
  (This reverses the earlier R5 "core owns Verbosity" decision, to follow new_c++.)
- **Mtilde** kept in `pairstat/StructureFactor` (new_c++ lacks it; FFH needs it).
  `GetKs` is in core, as in new_c++.
- **MultiStealthyPotential.{h,cpp}** brought into CCOptimization (new main needs it).
- **examples.cpp** re-adapted to the lean NumberVariance API
  (`MC_nv`→`MC_NV_Ensemble`, out-param `GetSamplingPoints`, +`<NumberVariance.h>`,
  `<RSA_Gen.h>`); it is scratch not built in new_c++ either.

Files only in `new_c++` and not needed by my builds were left out (e.g.
`Debye/`, `weighted/`, `PairStat/Percolation.h`, `FPT_Polygon/Spheres`,
`Trajectory`, `Potential/MonteCarlo.cpp`). All modules + consumers rebuild green.

## Module status

| Module | Tier | State | Builds here | Notes |
|--------|------|-------|-------------|-------|
| core | base | ✅ done | yes | smoke.out; owns `GetKs` (in PeriodicCellList) |
| pairstat | mid | ✅ done (lean) | yes | pairstat.out; `Mtilde` ported in from full SF |
| potential | mid | ✅ done | yes | objects only |
| models | mid | ✅ done | yes | GenerateConfigs from top-level |
| nonspherical | mid | ✅ done | yes | ParticleShapes + FormFactors + Dispersions (core-only) |
| inverse | high | ⛔ not started | — | needs nlopt (avail); (Gurobi opt) |
| tiling | high | 🟡 migrated, blocked | partial | sources split into `apps/Tessellation-based/` (see below); 3/5 objects compile, `Computation`/`main` blocked on lean-pairstat gaps + dead `LatticeGen` API |
| visualization | leaf | 🟡 impls migrated, 3/4 build | partial | headers `Plots.h`/`DisplaySpheres.h` + impls `math3d.{h,cpp}`, `glframe.h`, `Plots.cpp`, `DisplaySpheres.cpp`. `math3d.o`/`DisplaySpheres.o` compile (OpenGL `#ifdef`-guarded out); `Plots.o` **blocked** on misc `LatticeSumSolver` (not yet migrated). Drawing needs MathGL/OpenGL+freeglut (absent). |
| phonons | leaf | ✅ done | yes | `PhononFrequency.{h,cpp}` (+ `PlotPhononFrequencies.cpp`); both objects compile (core+potential+GSL; MathGL plotting guarded out) |
| lp | leaf | 🟡 migrated, blocked | yes (empty) | `TorquatoJiao.{cpp,h}` + vendored `gurobi_c*.h`; whole TU is `#ifdef USE_TJ` → empty object by default. A `USE_TJ=1` build needs glpk + Gurobi (absent). |
| misc | leaf | 🟡 headers only | n/a | `Discrepancy.h`, `ParameterSet.h`, `Interfaces.h` relocated into `modules/misc/`; header-only support shared by apps. Pending: the support/optimization library (`LatticeSumSolver`, `Optimize`, …) — see TODO. |

### Note on the lean-vs-full StructureFactor (resolved)
Two consumers (FFH, CCOptimization) needed `GetKs`/`Mtilde` that the lean
`pairstat` lacked. Resolution (R2 option b): `GetKs` already lives in core's
`PeriodicCellList`; only `Mtilde` was ported into `pairstat/StructureFactor.{h,cpp}`
(core-only, no nonspherical baggage). No consumer needed the dropped nonspherical
spectral-density functions.

## Consumer migration status

| Consumer | → | Deps | Builds here | State |
|----------|---|------|-------------|-------|
| examples (reference) | modules/examples | core, pairstat | yes | ✅ |
| CHI | modules/apps/CHI | core | yes | ✅ builds & runs |
| NR | modules/apps/NR | core, pairstat, models | yes | ✅ builds & runs (GetSamplingPoints adapted to lean API; now also builds examples.cpp) |
| FirstPassageTime | modules/apps/FirstPassageTime | core, nonspherical | yes | ✅ builds & links (runtime needs its dispersion data files) |
| CCOptimization | modules/apps/CCOptimization | core, pairstat, potential | yes | ✅ builds & runs (only main.cpp; heavy includes were dead #if 0) |
| FFH | modules/apps/FFH | core, pairstat (Mtilde), FFTW | yes | ✅ builds & links (FFTW_PREFIX → conda `pmp` env; runtime needs input data) |
| CPmanager (was EXC_formats) | modules/apps/CPmanager | core | yes | ✅ builds (ConfigPack ↔ txt manager; runtime needs args) |
| inverse (CC-inversion) | modules/apps/inverse | core, pairstat, potential, misc, visualization | objects compile | ✅ `main.o` + `CCInversionProblem.o` build. `main.cpp` is the entry point (the "CCInversion" command extracted from the legacy `Interfaces.h` cli(); defines `Verbosity` + `masterDir`). Linking a runnable binary still needs the visualization impls (MathGL/OpenGL). |
| Tessellation-based (tiling) | modules/apps/Tessellation-based | core, pairstat, models, misc, visualization | partial (3/5) | 🟡 **blocked**. `SpacefillingPacking.o`, `Voronoi.o`, `PolyPacking.o` compile; `Computation.o`/`main.o` do not (see split note below). `main.cpp` is the "Tiling" entry point extracted from the legacy `Interfaces.h` cli(). |

### Support-module split (visualization / misc)
The headers the CC-inversion app pulled in but that aren't part of any compute
module were moved out of `apps/inverse/` into two leaf header-only modules:
- `modules/visualization/` ← `Plots.h`, `DisplaySpheres.h` (the plan's `viz`
  tier, renamed). Declarations only; impls need MathGL/OpenGL.
- `modules/misc/` ← `Discrepancy.h`, `ParameterSet.h`, `Interfaces.h`.

Cross-includes (`Interfaces.h`→Plots/DisplaySpheres, `Plots.h`→ParameterSet) are
resolved by adding `-I$(misc) -I$(visualization)` to consumers. `config.mk` gained
`misc`/`visualization` path vars. `apps/inverse/main.cpp` was added and both app
objects still compile green.

### Tessellation-based split (tiling / Voronoi)
The legacy top-level tiling sources were migrated into `apps/Tessellation-based/`:
- `Voronoi.{h,cpp}` — relocated as-is (includes re-normalised to `<X.h>`; the
  unused `etc1.h` include dropped).
- `TilingProblem.cpp` (2100 lines) split into three units:
  - `SpacefillingPacking.{h,cpp}` — the `SpacefillingPacking` class +
    `getSpacefillingPacking` / `...2` drivers (legacy lines 16-340).
  - `Computation.{h,cpp}` — `ImLatGen` + the `Computation` task hierarchy
    (legacy lines 345-1362; classes are header-defined, all method bodies are
    inline, exactly as in the legacy single file).
  - `main.cpp` — `TilingProblemCLI` (legacy lines 1364-2100) + a small `main()`
    that is the **"Tiling" command** extracted from `misc/Interfaces.h` cli().
- `PolyPacking.{h,cpp}` — the `MakeMono2Poly*` family extracted from the retired
  `etc1.{h,cpp}` (only tiling consumes it, via `MakeMono2Poly_parallel`).

**Build status: blocked (3/5 objects).** `SpacefillingPacking.o`, `Voronoi.o`,
`PolyPacking.o` compile under default flags (the CGAL/voro++/MathGL branches are
`#ifdef`-guarded out — none of those libs are installed). `Computation.o` and
`main.o` do **not** compile, for reasons inherent to the legacy code, not the
split:
1. **lean-pairstat gaps (R2):** `Computation` needs `IsotropicSpectralDensity`
   /`IsotropicSpectralDensity_Lstorage` (spectral density) and `MC_nv` / `MC_vfv`
   (number- and volume-fraction-variance Monte-Carlo) that the lean `pairstat`
   dropped/renamed. Unblocking needs those ported into `pairstat`.
2. **dead `LatticeGen` API:** `ImLatGen` calls `LatticeGen(basis,int)`,
   `SetBasis`, `SetNumInEachSides`, `currPhi`, `rMax`, ... — all commented out in
   *every* surviving `GenerateConfigs.h` (top-level and `Models/`). The legacy
   `TilingProblem.cpp` was already stale against the current `GenerateConfigs`.
3. **heavy libs:** a fully-functional binary additionally needs CGAL + voro++ +
   MathGL (`make USE_CGAL=1 USE_VORO=1 USE_MATHGL=1`).

The migration (relocate + split + extract + entry point) is faithful and the
wiring is in place; resurrecting (1)/(2) is a separate porting effort.

### etc1.h retirement
Goal: remove the legacy grab-bag `etc1.{h,cpp}`. Usage audit — the only `etc1`
functions used by **compiled modular code** are `GetLattice_config` (CCInversion.h
+ CCInversionProblem.cpp) and `Find_LeastCommonDivider` (CCInversionProblem.cpp);
both have a single consumer (the CC-inversion app), so per the rule (infrequent →
move into the using source) they were relocated into
`apps/inverse/LatticeUtil.h`. The `MakeMono2Poly*` family is used by the tiling
app (`MakeMono2Poly_parallel` in `TilingProblemCLI`), so it was relocated into
`apps/Tessellation-based/PolyPacking.{h,cpp}` (see the Tessellation-based split).
All other `etc1` symbols (`GetLattice`, `ccConfig`,
`PHI`, `GetPatchedPoisson`, `WriteToContourf`, `ReadCoordinationNumber`,
`GetSimplePolyCrystals`, `iSoddyCircles`, `GetCantorSet`,
`ApplyCantorSet2Packing`, `LinealRemoval2D`,
`GetFibonacciChain`) are referenced only by uncompiled scratch
(`CCOptimization/DebugSimple.cpp`) or unmigrated legacy mains — none are needed by
the modular build. Result: **no `etc1.*` file exists in `modules/`**, and no
compiled TU includes `etc1.h`. (`Merge_Sort`/`MergeSort`, reached via etc1.h in the
legacy include chain, actually live in `Discrepancy.h`, which the app now includes
directly.)

## Unresolved / excluded source files (from the 5 migrated consumers)

Files present in a legacy consumer dir but **not** part of its migrated build,
with the reason. These are intentionally excluded (alternate mains / scratch /
dead), matching the legacy makefiles — they are not lost (legacy dirs untouched).

| Consumer | File | Reason unresolved |
|----------|------|-------------------|
| FFH | `FFH_padded.{h,cpp}` | alternate "padded" implementation; not in legacy build |
| FFH | `FFH_test.cpp` | alternate `main()` (test driver) |
| FFH | `sample1.cpp`, `sample2.cpp` | alternate `main()` scratch programs |
| FFH | `Dump.hpp` | helper header, not referenced by the built sources |
| NR | `distance-metrics.cpp` | scratch; `gaussian_pdf`/`compute_l2_and_kl` not declared in any NR header nor referenced by built sources |
| CCOptimization | `DebugSimple.cpp` | standalone scratch driver; not in legacy build |
| CHI | — | none (all sources resolved) |
| FirstPassageTime | — | none (all sources resolved) |

`nonspherical` module pulled only `ParticleShapes`, `FormFactors`, `Dispersions`
from the legacy `NonsphericalParticles/`. Still unresolved there (not needed by
the 5 consumers): `OrderMetrics.cpp`, `PeriodicNetworks.{h,cpp}`,
`SpectralDensity_Dispersion.{h,cpp}`, `Test_ParticleShapes.cpp`, `main.cpp`.

## Not yet migrated (legacy consumers remaining)
HSF, Statistics, QC, RandOrg, strong_contrast, SurfaceCoeff, THIRD_ORDER,
Example, NonsphericalParticles (as its own app), Models (CLI).
(Tessellation-based is now migrated under `apps/Tessellation-based/` — blocked,
see above.)

**Excluded from objectives (do not migrate):** Energy, Projections.

NR's `examples.cpp` is now built (added `RSA_Gen.h` to the `models` module with its
broken include guard fixed; adapted `MC_nv`→`MC_NV_Ensemble` and `GetSamplingPoints`
to the lean pairstat API).

## TODO — top-level `../c++/` sources not yet imported into `modules/`

Top-level `../c++/` sources, grouped by tier. Items migrated so far are ticked;
the rest are still only in the legacy tree (any `#include` of them in `modules/`
is inert — uncompiled scratch, behind `#ifdef VORONOI_INCLUDED`, or commented).

- **inverse**: ✅ `CCInversion.h`, `CCInversionProblem.{cpp,h}`
  migrated to `apps/inverse/` (with `LatticeUtil.h`). Still pending:
  `InvStatMechInterface.cpp`.
- **misc** (support / optimization library → `modules/misc/`): pending migration
  of `LatticeSumSolver.{cpp,h}`, `Optimize.{cpp,h}`, `ThreeBody.{cpp,h}`,
  `Fitting.{cpp,h}`, `Degeneracy.{cpp,h}`, `Solvers.h`,
  `withGaussianInvParameterization.h`, `Polynomial.h`. (These were previously
  grouped under "inverse library"; they belong in `misc`, which until now has
  held only relocated headers — so it would gain its first compiled sources.)
  ⚠ `LatticeSumSolver` here **blocks `visualization/Plots.o`** (the `Plots(...)`
  overloads take `std::vector<LatticeSumSolver*>` and dereference `->Terms`
  outside the MathGL guard).
- **misc headers (done)**: ✅ `Discrepancy.h`, `ParameterSet.h` → `misc/`;
  `Plots.h`, `DisplaySpheres.h` → `visualization/`; `Interfaces.h` → `misc/`.
  ✅ `etc1.{cpp,h}` retired (`GetLattice_config`/`Find_LeastCommonDivider` →
  `apps/inverse/LatticeUtil.h`; `MakeMono2Poly*` → `apps/Tessellation-based/
  PolyPacking`). Still pending: the `*.cpp` impl of `Discrepancy` (its header was
  relocated; `Plots.cpp`/`DisplaySpheres.cpp` are now migrated — see below).
- **tilings / Voronoi** → ✅ migrated to `apps/Tessellation-based/`
  (`TilingProblem.{cpp,h}` split + `Voronoi.{cpp,h}`; blocked, see above)
- **linear programming (Gurobi)** → ✅ migrated to `modules/lp/`
  (`TorquatoJiao.{cpp,h}` + vendored `gurobi_c*.h`; compiles empty, `USE_TJ=1`
  blocked on glpk/Gurobi)
- **phonons** → ✅ migrated to `modules/phonons/`
  (`PhononFrequency.{cpp,h}` + `PlotPhononFrequencies.cpp`; both objects build)
- **visualization** (impls) → ✅ migrated to `modules/visualization/`
  (`math3d.{h,cpp}`, `glframe.h`, `Plots.cpp`, `DisplaySpheres.cpp`; `Plots.o`
  blocked on misc `LatticeSumSolver`, the rest build)
- **CLIs / drivers**: `CollectiveCoordinateCLI.cpp`, `PairStatisticsCLI.cpp`
