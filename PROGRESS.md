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
`Debye/`, `PairStat/Percolation.h`, `FPT_Polygon/Spheres`,
`Trajectory`, `Potential/MonteCarlo.cpp`). All modules + consumers rebuild green.
(Note: `new_c++/Potential/MonteCarlo.cpp` was correctly left out — it is
ill-formed. The `MonteCarlo<>::Move/Anneal` bodies it was meant to hold were
instead restored **inline** into `modules/potential/MonteCarlo.h` from the
legacy `c++/Potential/` copy during the 2026-07-07 reorg; see the potential row
and `potential/CLAUDE.md`.)
(`new_c++/weighted/` has since been migrated as the `apps/weighted` consumer — see
its row in the consumer table below.)

## Module status

| Module | Tier | State | Builds here | Notes |
|--------|------|-------|-------------|-------|
| core | base | ✅ done | yes | smoke.out; owns `GetKs` (in PeriodicCellList) |
| pairstat | mid | ✅ done (lean) | yes | pairstat.out; `Mtilde` ported in from full SF. Doc comments reformatted 2026-06-30 to own-line block `/** \brief … */` + backslash tags (matches `apps/weighted`; `@p`→`\p`, `//<`→`///<`); 6 library objects rebuild clean, change verified comments-only |
| potential | mid | ✅ done | yes | **Full orchestrated code review 2026-07-22** (codex critic → Opus implementer → codex verifier, PASS): 29 findings fixed incl. physics bugs (LJ force σ-powers, MD initial acceleration, Nosé–Hoover half-velocity, AnalyzeData lag-0), RAII rebind in CellMove, covariant clone() overrides, ownership/copy-control, seed-scheme, latent-template fixes gated by a forced-instantiation TU; all "(CLAUDE)" markers finalized (14 paper-eq citations left for author eyeball — see `potential/CLAUDE.md`); documented param names restored. objects only. Header reorg 2026-07-07 (see `potential/CLAUDE.md`): cooling schedules split into `CoolingSchedule.h`; CC potentials split into base + `CCPotential_varBox.h` + `CCPotential_Kernels.h`; `Move.h`→CC vestigial include dropped; missing `MonteCarlo<>::Move/Anneal` template bodies restored from legacy (were declared-but-undefined → latent link error). **Warning-cleanup 2026-07-08: `make objs` now 0 warnings** (also strict `-Wall -Wextra -Wpedantic`) — fixed the C++20-ill-formed `ParallelTempering` ctor template-id, 2 `"\%"` iostream escapes in `StructureOptimization.cpp`, and ~30 `-Wunused-parameter` sites (param names commented out); comments-only + rebuild verified |
| models | mid | ✅ done | yes | GenerateConfigs from top-level. **Full code review 2026-07-22** (orchestrated: codex critic → Opus implementer → codex verifier, final PASS): all 16 findings from the 2026-07-06 audit fixed + ~28 new ones — RandomOrganization min-image/active-test/OMP-race fixes, `RSA_Gen.h` ODR `inline` fix, `unique_ptr` ownership (HSF, GenerateConfigCLI), RNG reproducibility (deterministic per-thread seeding), `PhiC`/`Rc` consolidated (1D = 0.7476 Rényi), documented param names restored per Doxygen rule. Accepted-as-is by decision: R07 partial (no model dispatch fabricated), N16 (numPrts-based rescaling documented). Both bars 0 warnings; consumers NR/examples/Tessellation rebuild green. README gained model definitions + citations (Zhang–Torquato PRE 88, 053312 (2013); Tjhung–Berthier PRL 114, 148301 (2015)). Details: `models/CLAUDE.md`. **Merged with legacy `c++/Models/` + reorganized into per-model subdirs 2026-07-22** (orchestrated, final verifier PASS): `RSA/`, `HSF/` (+ ported `HSF_Gen.h` — Tessellation `main.cpp` blocker now gone), `Pert/GenPerturb.h` + `Lattices/Lattice.{h,cpp}` (split out of the GenerateConfigs monolith), `RandOrg/`; flat includes preserved via `config.mk` `models_subdirs_INC` (zero consumer churn, `git mv` renames staged). Ported: `ConfigGen` param/stream ctors, `Echo<T>`, `HSF_Gen`. Pruned dead `LoadConfig`/`GenerateConfigCLI` (N22/R07 → resolved by the working `apps/Models` CLI). All 53 "(CLAUDE)" doc-comment placeholders audited by codex and finalized (15 corrected, markers removed). **Both passes committed as `abb1fbe`** (21 files, renames preserve history). |
| nonspherical | mid | ✅ done | yes | ParticleShapes + FormFactors + Dispersions (core-only). **+ `SpectralDensity_Dispersion.{h,cpp}` (2026-07-16)**: orientationally-averaged χ_V(k) of nonspherical dispersions (the nonspherical analog of pairstat's spherical `SpectralDensity`; relocated from legacy `NonsphericalParticles/`, `GetKs` now from core so no pairstat dep). Fixed a missing-`return` UB in it. See `nonspherical/CLAUDE.md`. |
| inverse | high | ➡️ built as `apps/inverse` | yes (app) | no standalone `modules/inverse/` library — the CC-inversion work is the **`apps/inverse`** consumer (builds & runs, see below). needs nlopt (avail); (Gurobi opt) |
| tiling | high | ✅ compiles | yes (objs) | sources split into `apps/Tessellation-based/` (see below); **all 5 objects compile** (2026-06-23: `ImLatGen` rewritten to the current `LatticeGen` API; spectral-density + `MC_vfv` ported into pairstat). A *runnable* binary still needs CGAL/voro++/MathGL. |
| visualization | leaf | ✅ all objects build | yes | headers `Plots.h`/`DisplaySpheres.h` + impls `math3d.{h,cpp}`, `glframe.h`, `Plots.cpp`, `DisplaySpheres.cpp`. **Plots.o unblocked 2026-07-23** (final audit X3-01: the two `Plots()` definitions gained the same `LATTICESUMSOLVER_INCLUDED` guard as their declarations) — default `make` now compiles all objects, 0 warnings. Drawing still needs MathGL/OpenGL+freeglut (absent); those paths + the LatticeSumSolver overloads remain static-review-only. 14 audit bug fixes 2026-07-23 (see `visualization/CLAUDE.md`). |
| phonons | leaf | ✅ done | yes | `PhononFrequency.{h,cpp}` (+ `PlotPhononFrequencies.cpp`); both objects compile (core+potential+GSL; MathGL plotting guarded out) |
| lp | leaf | 🟡 migrated, blocked | yes (empty) | `TorquatoJiao.{cpp,h}` + vendored `gurobi_c*.h`; whole TU is `#ifdef USE_TJ` → empty object by default. A `USE_TJ=1` build needs glpk + Gurobi (absent). |
| misc | leaf | 🟡 headers | n/a | `Discrepancy.h`, `ParameterSet.h`, `Interfaces.h`, **`LatticeUtil.h`** (promoted from apps/inverse, 2026-06-29) and **`misc.h`** (`SimpleHarmonic`/`SimpleHarmonicF`/`HS_pot`, moved out of apps/inverse 2026-06-29) in `modules/misc/`; header-only support shared by apps. ⚠ `misc.h` makes misc depend on **potential** (the three classes are `Potential` subclasses). Pending: the support/optimization library (`LatticeSumSolver`, `Optimize`, …) — see TODO. |

### Note on the lean-vs-full StructureFactor (resolved)
Two consumers (FFH, CCOptimization) needed `GetKs`/`Mtilde` that the lean
`pairstat` lacked. Resolution (R2 option b): `GetKs` already lives in core's
`PeriodicCellList`; only `Mtilde` was ported into `pairstat/StructureFactor.{h,cpp}`
(core-only, no nonspherical baggage). No consumer needed the dropped *nonspherical*
spectral-density functions.

Update (2026-06-23): the tiling app needed the **spherical** spectral density
of sphere packings, so `SpectralDensity` / `SpectralDensity_Lstorage` and the
ensemble `IsotropicSpectralDensity` / `IsotropicSpectralDensity_Lstorage` were
ported into `pairstat/StructureFactor.{h,cpp}` (spherical only — still no
nonspherical/USE_SHAPES baggage), and the local volume-fraction variance
`MC_vfv` (+ file-local `I_1`/`I_config`) into `pairstat/NumberVariance.{h,cpp}`.

## Final audit (2026-07-23)

Whole-tree cross-family audit (4 codex shards + 3 Opus shards, 198 findings)
followed by an Option-A fix round (4 Opus implementers, opposite-family
confirmation): **126 bug-band fixes**, design/style findings documented in each
unit's `CLAUDE.md` (`## 2026-07-23 final audit` sections; files created for
core, pairstat, misc, NR, PairStatCLI, CPmanager, inverse, weighted).
Highlights: core `GetNeighborList` cache corruption; pairstat `vSkBin`-under-
Verbosity, ensemble-g₂ normalization, H_v sampling, Haar rotations, lineal-path
tail; `misc/ParameterSet.h` copy/move/delete[] criticals; TJ jammer lifecycle;
**FFH zero-iteration regression (C1-01) fixed** — default `Run()` converges
again; `NetworkBV.out` live B_V path (bias/OOB/format); NR L2-norm
underestimate; visualization Plots.o unblocked. **Numerics-affecting fixes are
tabled per-area in `_multiagent/tasks/final-audit/workers/fix/resultF1-F4.md`**
— outputs of g₂/S(k)/H_v/vfv/L2/B_V/FFH-ε analyses can legitimately change.
Author call still pending: X1-34 (`ClusterSize` radius-vs-diameter doc).
Integration build after fixes: all 9 modules + 13 apps + examples, 0 warnings.
README.md gained the per-app purpose table.

## Consumer migration status

| Consumer | → | Deps | Builds here | State |
|----------|---|------|-------------|-------|
| examples (reference) | modules/examples | core, pairstat | yes | ✅ |
| CHI | modules/apps/CHI | core | yes | ✅ builds & runs |
| NR | modules/apps/NR | core, pairstat, models | yes | ✅ builds & runs (GetSamplingPoints adapted to lean API; now also builds examples.cpp) |
| FirstPassageTime | modules/apps/FirstPassageTime | core, nonspherical | yes | ✅ builds & links (runtime needs its dispersion data files) |
| CCOptimization | modules/apps/CCOptimization | core, pairstat, potential | yes | ✅ builds & runs (only main.cpp + MultiStealthyPotential.cpp; heavy includes were dead #if 0). Re-verified 2026-07-08 against the reorganized `potential` headers. **Full code review 2026-07-22** (orchestrated, PASS): 16 findings fixed — RepulsiveCC AllForce double-count, MultiStealthy species bookkeeping (indices/duplicates/radii/coverage validation), end()-deref, format-string UB, N_load underflow now a loud error, unique_ptr drivers; "(CLAUDE)" markers finalized (DebugSimple.cpp excluded). Details: `apps/CCOptimization/CLAUDE.md`. Relinks `soft_core_stealthy.out`, 0 warnings both bars. |
| FFH | modules/apps/FFH | core, pairstat (Mtilde), FFTW | yes | ✅ builds & links (FFTW_PREFIX → conda `pmp` env; runtime needs input data) |
| CPmanager (was EXC_formats) | modules/apps/CPmanager | core | yes | ✅ builds (ConfigPack ↔ txt manager; runtime needs args) |
| inverse (CC-inversion) | modules/apps/inverse | core, pairstat, potential, misc | ✅ builds & runs | ✅ **COMPLETE: `make` → runnable `CCInversion.out`; example runs; 0 compiler warnings in its own sources** (verified 2026-06-29; re-verified 2026-07-08 against the reorganized `potential` headers — relinks `CCInversion.out`, still 0 warnings). **No MathGL/OpenGL needed** — the "CCInversion" command doesn't touch the visualization impls (`-I$(visualization)` is on the include path but no viz objects are linked). `main.cpp` is the entry point (extracted from the legacy `Interfaces.h` cli(); defines `Verbosity` + `masterDir`). Has `README.md` + `run_example.sh`. |
| Tessellation-based (tiling) | modules/apps/Tessellation-based | core, pairstat, models, visualization | ✅ **runnable** (`USE_VORO`) | ✅ **`make USE_VORO=1` → `TessellationCLI.out`, 0 warnings / 0 errors** (2026-07-16). All five objects build; see the split note below for the fixes. `TessellationAnalysisCLI.cpp` is the default (config-reading) CLI; `main.cpp` (the generator-coupled "Tiling" entry point) still needs the unmigrated `models/HSF_Gen.h`. **Decoupled from `misc`** (2026-07-16): the sole `misc/Discrepancy` symbol used (`Merge_Sort`) was copied file-locally into `Voronoi.cpp` (anonymous namespace), so `-I$(misc)` / `Discrepancy.o` are gone; `misc/Discrepancy.*` untouched (still serve `apps/inverse`). `README.md` now also documents **how to run** the analysis-only CLI: it is `stdin`-driven (no analysis argv; `argv[1]` = working dir only), with a concrete `TessellationCLI.out` command-script example + a command-tree / task-parameter reference; a "CLI drivers" table makes explicit that `TessellationCLI.out` **is** the config-reading analysis-only CLI (vs. the un-built generator-coupled `main.cpp`). |
| weighted (marked pt stats) | modules/apps/weighted | core | ✅ builds & runs | ✅ **COMPLETE (2026-06-30).** Builds & runs `weighted.out`, warning-clean (depends on `core` only + header-only Boost/OpenMP). A simplified `PairStatisticsCLI` for **weighted/marked** point patterns: weighted structure factor (`SkComputation`), weighted two-/single-sided `g2` (`g2Computation`), and local variance of the weights (`LocalVariance`). Ships its own copies of `Computation`/`PairCorrelation`/`LocalVariance`/`WeightedStructureFactor` (distinct from the `pairstat` module's same-named files; they operate on `PeriodicCellList<std::complex<double>>`). Source = `new_c++/weighted/`. Used for the PRX 16, 011042 (2026) "Hyperuniformity of Weighted Particle Systems" paper. Doxygen comments added 2026-06-29. `README.md` (re)created 2026-06-30: documents the PRX citation, that the command language mirrors `apps/PairStatCLI`, and the weighted input format — `ReadConfigPack <configpack-prefix> <weights-prefix(printf %04d)> <column>` over `MarkedPointPattern = ConfigurationPack (.ConfigPack/.ConfigPackIndex) + per-config weights txt`. `test.sh`/`test2.sh` revised into runnable example scripts that mirror the README and reproduce the bundled `test1_*`/`test2_*`/`test3_*` outputs (Sk on 3D MRJ; LocalVariance ±`UseParticleCenters` and two-/single-sided g2 on 2D RSA); both pass `bash -n`. Doc-comment spacing normalized 2026-06-30 (a blank line now precedes each function's leading comment in `Computation.h`/`Computation.cpp`/`main.cpp`); whitespace-only, rebuilds green. Compiles **warning-clean** 2026-06-30 (`make` from a clean tree → 0 warnings/0 errors): removed all ~35 `-Wsign-compare` + `-Wunused-variable` + lambda-capture warnings across the 5 sources — signed loop counters widened to `size_t`, 4 dead `mean`/`R` locals dropped, and the global `Rmin_local_variance` no longer captured by a lambda. Doc comments unified 2026-06-30 to block `/** \brief … */` style with backslash Doxygen tags (`\brief`/`\param`/`\return`/`\p`) across all 9 sources (`///`-runs folded into blocks; at-sign tags converted); rebuilds warning-clean. |

| CollectiveCoordinateCLI | modules/apps/CollectiveCoordinateCLI | core, pairstat, potential, models (+ misc/visualization headers) | ✅ builds & runs | ✅ **COMPLETE (2026-07-22, orchestrated, verifier PASS).** `CollectiveCoordinateCLI.out` — **Ge Zhang's** stealthy-hyperuniform ground-state CLI (used for Sci. Rep. 6, 36963 (2016); PRE 92, 022120 (2015); PRE 92, 022119 (2015); PRX 5, 021020 (2015) — see its README). Faithful adapted copy of legacy top-level `CollectiveCoordinateCLI.cpp` (2694→2314 lines, CRLF preserved): `Debug` command dropped (sole consumer of unmigrated Solvers/Degeneracy), `::Plot`/`PlotDirectionalStructureFactor2D` = app-local no-op stubs, `PlotFunction_Grace` = byte-identical legacy copy (WLMC `.agr` outputs work), `USE_EIGEN` off (built-in stubs). main.cpp replicates legacy startup (Verbosity, ProgramStart, **TimeLimit** — 0 would exit(0) in RelaxStructure_NLOPT — nlopt_srand, omp nesting, precision 10). Inherited-defect fixes: null-state guards on dependent commands (AutoRescale correctly pConfig-only), ReadProblemFile clean null state, temp-and-swap generation. `run_example.sh` end-to-end: 2D N=40 χ=0.3 → E ~1.7e-21. 0 warnings both bars; zero module edits. |
| Models (CLI) | modules/apps/Models | core, models | ✅ builds & runs | ✅ **COMPLETE (2026-07-22).** `models.out` — the model-generation CLI adapted from legacy `c++/Models/CLI.cpp` (which was the *working* implementation of what the pruned `GenerateConfigCLI`/`LoadConfig` placeholders promised). Models: Poisson, RSA, HSF, Perturbations; output ConfigPack or per-file SpherePacking; `ReadConfigPack` input for Perturbations. **stdin token order = legacy protocol** (verified against the legacy README/stream ctors; end-to-end RSA/HSF/Poisson/Perturbations script runs tested). Legacy CLI bugs fixed in adaptation: uninitialized model pointer on unknown names, prefix-used-as-format-string, unvalidated interactive input (d ≤ MaxDimension, counts ≤ 1e8, finite ρ), Perturbations+SpherePacking combo now rejected, `unique_ptr` ownership, `size_t` loops. Has `README.md` + `run_example.sh` (both corrected 2026-07-22 after a user-caught stale token order — README/example predated the legacy-protocol remediation of `main.cpp`). RSA at `Verbosity<=5` logs insertion progress to a lazily-created `DetailedLog.txt` (legacy behavior, kept by decision). Committed in `abb1fbe`. 0 warnings both bars. |
| Statistics (n-body $g_n$) | modules/apps/Statistics | core (gn.out); core, models (tests/gen.out) | ✅ builds & runs | ✅ **COMPLETE (2026-10-01).** Migrated from `new_c++/Statistics/` (gn, tau_n, main; the test-input generator `gen.cpp` and `check_url_gn.py` live in `tests/` with their own makefile). Doxygen on every function; 0 warnings both bars; ASan/UBSan clean. `gn` kernel rewritten for speed (sorted window scan, one parallel region per config over owned x_21 slabs, touched-bin accumulation; thread-count independent): 40–2500× on the sanity cases, ~6× per config for n=4 at Rmax=10. Inherited defects fixed: half-filled outermost bins, rho=1 normalization, non-integer Rmax/dx abort, `gn_1D … n` flag ignored, `Compute_tau_n_excluded` result[-1]. Validated by `tests/check_url_gn.py` against the exact g_n of the 1D cloaked URL (a=1,2; n=2,3,4: PASS); migrated interior bins identical to legacy. Details: `apps/Statistics/README.md`. |
| OrderMetric_PeriodicNetwork | modules/apps/OrderMetric_PeriodicNetwork | core, nonspherical (+ pairstat, models for the network target) | ✅ builds & links | ✅ **COMPLETE (2026-07-16).** Migrated from legacy `NonsphericalParticles/`. **Two runnable targets** (the folder bundles two non-interacting programs): `GetChivDispersions.out` ← `main.cpp` — the orientationally-averaged spectral density χ_V(k) of nonspherical **dispersions** (`IsotropicSpectralDensity_Dispersion` + form-factor tests; **NOT** an order-metric driver despite the folder name), links core + nonspherical only; and `NetworkBV.out` ← a new thin `bv_main.cpp` wrapping `PeriodicNetworks`'s `ComputeBV_CLI()` (surface-area coefficient B̄_V / volume-fraction variance of periodic two-phase media), links core + nonspherical + pairstat + models. `SpectralDensity_Dispersion.{h,cpp}` was **relocated to the `nonspherical` module** (see its row) and the two legacy `SpectralDensity_nonspherical` calls in `PeriodicNetworks.cpp` rewritten to it. Needs `-DUSE_SHAPES`. Behavior-preserving c++17 fixes: a GNU default-arg-temporary-bound-to-ref hard error, ~30 `-Wsign-compare`, unused-vars, maybe-uninitialized. 0 warnings / 0 errors both targets. `PeriodicNetworks`'s historical order-metric `main()` was `Test_ParticleShapes.cpp` (**not** migrated). Runtime: `GetChivDispersions.out` needs `.Dispersion` input files; `NetworkBV.out` is interactive stdin. Has `CLAUDE.md`. |

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
- `Voronoi.{h,cpp}` — relocated as-is (includes re-normalized to `<X.h>`; the
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

**Build status: all 5 objects compile (2026-06-23).** Previously `Computation.o`
and `main.o` were blocked; both now build. What was fixed:
1. **lean-pairstat gaps (R2) — RESOLVED by porting into `pairstat`.**
   `IsotropicSpectralDensity` / `IsotropicSpectralDensity_Lstorage` (+ the
   single-config `SpectralDensity` / `SpectralDensity_Lstorage`) were ported into
   `pairstat/StructureFactor.{h,cpp}`, and `MC_vfv` (+ the `I_1`/`I_config`
   helpers) into `pairstat/NumberVariance.{h,cpp}` — same approach used earlier
   for `Mtilde`. `MC_nv` was switched to the lean equivalent `MC_NV_Ensemble`.
2. **dead `LatticeGen` API — RESOLVED by rewriting `ImLatGen`.** It no longer
   uses `LatticeGen(basis,int)`/`SetBasis`/`SetNumInEachSides`/`currPhi`/`N`/`rMax`;
   it builds a 1-particle cubic unit cell, keeps its own tile-count + radius,
   tiles, and deletes random vacancies (Computation.h).
3. **missing generators — RESOLVED.** `HSF_Gen` migrated from `c++/HSF/HSF_Gen.h`
   into `modules/models/HSF_Gen.h`; `PoissonPt` (absent from the modular
   `GenerateConfigs`) reintroduced as a small generator in `main.cpp`; `main.cpp`
   gained `#include <RSA_Gen.h>` (for `RSA_Gen`/`PhiC`).
4. **heavy libs (still required for a runnable binary):** the CGAL/voro++/MathGL
   branches stay `#ifdef`-guarded out; a fully-functional binary needs
   `make USE_CGAL=1 USE_VORO=1 USE_MATHGL=1` once those libs are installed.

The migration is faithful and all objects compile under the default flags;
the regression set (core/pairstat/potential/models + examples/NR/PairStatCLI/
CCOptimization) still builds and links green after the pairstat additions.

### etc1.h retirement
Goal: remove the legacy grab-bag `etc1.{h,cpp}`. Usage audit — the only `etc1`
functions used by **compiled modular code** are `GetLattice_config` and
`Find_LeastCommonDivider` (CCInversionProblem.cpp). They were first relocated into
`apps/inverse/LatticeUtil.h` (single consumer), then **promoted to
`modules/misc/LatticeUtil.h` on 2026-06-29** once a second consumer appeared
(`misc.h`'s `SimpleHarmonicF`) — consumers now `#include <LatticeUtil.h>` via
`-I$(misc)`. The `MakeMono2Poly*` family is used by the tiling
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
legacy include chain, lived in `misc/Discrepancy`; the app used only `Merge_Sort`,
so on 2026-07-16 it was **copied file-locally into `Voronoi.cpp`** and the app was
fully decoupled from `misc`.)

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

`nonspherical` module pulled `ParticleShapes`, `FormFactors`, `Dispersions`
(and, 2026-07-16, `SpectralDensity_Dispersion.{h,cpp}`) from the legacy
`NonsphericalParticles/`. Migrated 2026-07-16 into `apps/OrderMetric_PeriodicNetwork/`:
`PeriodicNetworks.{h,cpp}` + `main.cpp` (see that consumer row). Still unresolved
there: `OrderMetrics.cpp` (a 3-line stub, not built) and `Test_ParticleShapes.cpp`
(the historical order-metric `main()`; 3200-line test scaffold, not migrated).

## Not yet migrated (legacy consumers remaining)
HSF, Statistics, QC, RandOrg, strong_contrast, SurfaceCoeff, THIRD_ORDER,
Example.
(Models (CLI) is now migrated as `apps/Models` — see the consumer table. The
`c++/Models/` per-model subdir sources are merged into `modules/models/`.)
**Follow-up unblocked:** `apps/Tessellation-based/main.cpp` (generator-coupled
"Tiling" CLI) needed `HSF_Gen.h`, which now exists in `models/HSF/` — wiring it
into the Tessellation makefile is a pending follow-up (excluded from the
2026-07-22 merge task by decision).
(Tessellation-based is now migrated under `apps/Tessellation-based/` — compiles
(all 5 objects), see above.)
(NonsphericalParticles: its χ_V driver `main.cpp` + `PeriodicNetworks.{h,cpp}` are
now migrated under `apps/OrderMetric_PeriodicNetwork/` (2026-07-16); only the
`Test_ParticleShapes.cpp` order-metric driver remains legacy-only.)

**Excluded from objectives (do not migrate):** Energy, Projections.

NR's `examples.cpp` is now built (added `RSA_Gen.h` to the `models` module with its
broken include guard fixed; adapted `MC_nv`→`MC_NV_Ensemble` and `GetSamplingPoints`
to the lean pairstat API).

## TODO — top-level `../c++/` sources not yet imported into `modules/`

Top-level `../c++/` sources, grouped by tier. Items migrated so far are ticked;
the rest are still only in the legacy tree (any `#include` of them in `modules/`
is inert — uncompiled scratch, behind `#ifdef VORONOI_INCLUDED`, or commented).

- **inverse**: ✅ **complete.** `CCInversion.h`, `CCInversionProblem.{cpp,h}`
  migrated to `apps/inverse/`; builds & runs `CCInversion.out`. `CCInversion.h`
  now holds only the `FrankInvertCCPotential*` family —
  `SimpleHarmonic`/`SimpleHarmonicF`/`HS_pot` moved to `misc/misc.h` and
  `LatticeUtil.h` promoted to `misc/` (2026-06-29).
  (`InvStatMechInterface.cpp` is **not** part of the CC-inversion app — different
  subject — so it is not tracked here.)
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
  `misc/LatticeUtil.h` (2026-06-29; first lived in `apps/inverse/`);
  `MakeMono2Poly*` → `apps/Tessellation-based/PolyPacking`). ✅ `misc/misc.h`
  added (the `SimpleHarmonic`/`SimpleHarmonicF`/`HS_pot` potentials moved out of
  `apps/inverse/CCInversion.h`). Still pending: the `*.cpp` impl of `Discrepancy`
  (its header was relocated; `Plots.cpp`/`DisplaySpheres.cpp` are now migrated).
- **tilings / Voronoi** → ✅ migrated to `apps/Tessellation-based/`
  (`TilingProblem.{cpp,h}` split + `Voronoi.{cpp,h}`; compiles, see above)
- **linear programming (Gurobi)** → ✅ migrated to `modules/lp/`
  (`TorquatoJiao.{cpp,h}` + vendored `gurobi_c*.h`; compiles empty, `USE_TJ=1`
  blocked on glpk/Gurobi)
- **phonons** → ✅ migrated to `modules/phonons/`
  (`PhononFrequency.{cpp,h}` + `PlotPhononFrequencies.cpp`; both objects build)
- **visualization** (impls) → ✅ migrated to `modules/visualization/`
  (`math3d.{h,cpp}`, `glframe.h`, `Plots.cpp`, `DisplaySpheres.cpp`; `Plots.o`
  blocked on misc `LatticeSumSolver`, the rest build)
- **CLIs / drivers**: ✅ `CollectiveCoordinateCLI.cpp` → migrated as
  `apps/CollectiveCoordinateCLI` (2026-07-22; see its consumer row). Still
  pending: `PairStatisticsCLI.cpp`
