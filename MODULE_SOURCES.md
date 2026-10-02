# Module source provenance

This file records, per file present in `modules/`, **which legacy copy was
chosen as canonical** and why. It is regenerated from the *current* folder/file
layout (not just the originally-planned one) and cross-checked by content:
every source file under `modules/` was MD5-hashed and line-similarity-scored
against every same-named file in the two legacy trees, so the "From" column is
the measured closest source, not a guess.

Policy: *prefer the leaner subdirectory "module" copies* over the legacy flat
top-level copies, **except** where `new_c++/` provides an authoritative newer
version (then that wins) or where the subdir copy is a stub / missing features
(then the top-level copy wins). Every exception is flagged.

Legend: ✅ chosen per policy · ⚠️ exception to the default policy · 🆕 no legacy
source (authored for the modular build, or an extraction/split of legacy code).

## Legacy source roots

Two parallel legacy trees were used as sources (both left untouched, per R1):

| Root | What it is |
|------|------------|
| `../c++/` | The legacy codebase: a **flat top-level** copy of every file *plus* per-consumer subdirs (`cores/`, `PairStat/`, `Potential/`, `Models/`, `NonsphericalParticles/`, `CHI/`, `FFH/`, …). |
| `../new_c++/` | A **newer authoritative** parallel tree (subset of consumers: `cores/`, `PairStat/`, `Potential/`, `CCOptimization/`, `NR/`, `FirstPassageTime/`, `CPmanager/`, `weighted/`, …). Where a file exists here it is preferred (the "Synced to new_c++" pass in `PROGRESS.md`). |

Most `modules/` files are **not byte-identical** to their source: includes were
re-normalized `#include "../X.h"` → `#include <X.h>` (rule R4), so a file is
"from" the source it is *closest* to after that rewrite.

---

# Library modules

## core/
Source: **`new_c++/cores/`** ✅ (authoritative newer sync; closest match for all 8 files).

| File | From | Notes |
|------|------|-------|
| GeometryVector.{h,cpp} | new_c++/cores/ | include-normalized |
| PeriodicCellList.{h,cpp} | new_c++/cores/ | owns `GetKs` (used by FFH/CCOptimization) |
| RandomGenerator.{h,cpp} | new_c++/cores/ | |
| etc.{h,cpp} | new_c++/cores/ | globals (`Verbosity` etc.); each `main` defines `Verbosity` (new_c++ convention) |
| smoke.cpp | 🆕 | authored smoke-test driver (no legacy original) |

## pairstat/
Source: **`new_c++/PairStat/`** ✅ (lean, decoupled; newer than the `c++/PairStat/` and flat copies).

| File | From | Notes |
|------|------|-------|
| PairCorrelation.{h,cpp} | new_c++/PairStat/ | `.cpp` byte-identical to source |
| StructureFactor.{h,cpp} | new_c++/PairStat/ | ⚠️ `Mtilde` retained from the full legacy SF (new_c++ lacks it; FFH needs it). Still drops the *nonspherical* extras. |
| NumberVariance.{h,cpp} | new_c++/PairStat/ | lean API (`MC_NV_Ensemble`, out-param `GetSamplingPoints`) |
| ContactStatistics.{h,cpp} | new_c++/PairStat/ | `.cpp` byte-identical |
| SpectralDensity.{h,cpp} | c++/ (flat, spherical only) | 🆕 spherical spectral density `SpectralDensity[_Lstorage]` / `IsotropicSpectralDensity[_Lstorage]` ported from the full legacy StructureFactor (no nonspherical/USE_SHAPES baggage). Ported 2026-06-23 for the tiling app; subsequently split into this dedicated TU. |
| VolumeFractionVariance.{h,cpp} | c++/ (flat) | 🆕 local volume-fraction variance `MC_vfv` (+ file-local `I_1`/`I_config`) ported from the full legacy NumberVariance. Ported 2026-06-23 for the tiling app; subsequently split into this dedicated TU. |

Note: the legacy `PairStat/Computation.{h,cpp}` is **not** part of this module —
it migrated into `apps/PairStatCLI/` (the CLI consumer) instead.

## potential/
Source: **`new_c++/Potential/`** ✅, with two flat-only files from `c++/`.

| File | From | Notes |
|------|------|-------|
| Potential.{h,cpp} | new_c++/Potential/ | `.cpp` byte-identical |
| MD_System.{h,cpp} | new_c++/Potential/ | |
| StructureOptimization.{h,cpp} | new_c++/Potential/ | |
| MonteCarlo.h, Move.h, CollectiveCoordinatePotential.h, AnalyzeData.h | new_c++/Potential/ | header-only |
| LJPotential.{h,cpp} | **c++/** (flat) | ⚠️ no subdir / new_c++ copy; `.cpp` byte-identical |
| MC_System.{h,cpp} | **c++/** (flat) | ⚠️ no subdir / new_c++ copy |

## models/
Default policy would pick `Models/`, but the flat top-level copies win throughout:

| File | From | Notes |
|------|------|-------|
| GenerateConfigs.{h,cpp} | **c++/** (flat) | ⚠️ EXCEPTION: `Models/` copy is a stub (similarity 0.96 vs flat, 0.20 vs `Models/`) |
| HardSphereFluids.{h,cpp} | c++/ (flat) | `.cpp` byte-identical; no clean subdir module |
| RandomOrganization.{h,cpp} | c++/ (flat) | `.cpp` byte-identical; closer than `RandOrg/` |
| RandomSequentialAddition.{h,cpp} | c++/ (flat) | tied with `Models/RSA/`; flat chosen |
| RSA_Gen.h | c++/ (flat) | broken include guard fixed during migration |
| HSF_Gen.h | c++/HSF/ | 🆕 migrated 2026-06-23 (header-only `ConfigGen` wrapping the HSF sampler); includes normalized, methods made `inline`. Needed by the tiling app. |

## nonspherical/
Source: **`c++/NonsphericalParticles/`** ✅ (only location).

| File | From | Notes |
|------|------|-------|
| ParticleShapes.{h,cpp} | c++/NonsphericalParticles/ | |
| FormFactors.cpp | c++/NonsphericalParticles/ | |
| Dispersions.{h,cpp} | c++/NonsphericalParticles/ | |

Pulled only these three; `OrderMetrics`, `PeriodicNetworks`,
`SpectralDensity_Dispersion`, `Test_ParticleShapes`, `main.cpp` left in legacy.

## phonons/
Source: **`c++/`** (flat) ✅ (only location).

| File | From | Notes |
|------|------|-------|
| PhononFrequency.{h,cpp} | c++/ (flat) | `.cpp` byte-identical; MathGL plotting `#ifdef`-guarded |
| PlotPhononFrequencies.cpp | c++/ (flat) | |

## lp/
Source: **`c++/`** (flat) ✅. Whole TU is `#ifdef USE_TJ` (Gurobi/glpk absent → empty object).

| File | From | Notes |
|------|------|-------|
| TorquatoJiao.{h,cpp} | c++/ (flat) | linear-programming jamming |
| gurobi_c.h, gurobi_c++.h | c++/ (flat) | vendored Gurobi headers, byte-identical |

## misc/
Source: **`c++/`** (flat) ✅. Header-only support shared by apps.

| File | From | Notes |
|------|------|-------|
| Discrepancy.h | c++/ (flat) | also provides `Merge_Sort`/`MergeSort` |
| ParameterSet.h | c++/ (flat) | |
| Interfaces.h | c++/ (flat) | legacy `cli()` dispatcher (apps extract their entry points from it) |
| LatticeUtil.h | c++/etc1.h (via apps/inverse) | 🆕 `GetLattice_config` + `Find_LeastCommonDivider`; first in `apps/inverse/`, **promoted here 2026-06-29** (2nd consumer: `misc.h`). core-only. |
| misc.h | apps/inverse/CCInversion.h (orig. c++/CCInversion.h) | 🆕 `SimpleHarmonic`/`SimpleHarmonicF`/`HS_pot` potentials moved out of the inverse app 2026-06-29 (generic `Potential` subclasses; weren't used by any compiled source). ⚠ makes misc depend on **potential**; includes `<LatticeUtil.h>`. |

Pending (not yet migrated): the `.cpp` impls and the support/optimization library
(`LatticeSumSolver`, `Optimize`, `ThreeBody`, `Fitting`, `Degeneracy`, …).
`LatticeSumSolver` is what currently **blocks `visualization/Plots.o`**.

## visualization/
Source: **`c++/`** (flat) ✅. Drawing impls need MathGL/OpenGL (guarded out by default).

| File | From | Notes |
|------|------|-------|
| Plots.{h,cpp} | c++/ (flat) | `Plots.o` blocked on misc `LatticeSumSolver` |
| DisplaySpheres.{h,cpp} | c++/ (flat) | compiles (OpenGL `#ifdef`-guarded) |
| math3d.{h,cpp} | c++/ (flat) | compiles |
| glframe.h | c++/ (flat) | header-only |

## examples/
| File | From | Notes |
|------|------|-------|
| g2_sk_demo.cpp | 🆕 | authored reference consumer (g(r)/S(k) demo); no legacy original |

---

# Application drivers — `apps/`

Each app is a migrated legacy consumer (its own sources copied in, includes
re-normalized). "From" is the legacy consumer dir the sources came from.

## apps/CCOptimization/
Source: **`new_c++/CCOptimization/`** ✅ (main + potentials).

| File | From | Notes |
|------|------|-------|
| main.cpp | new_c++/CCOptimization/ | similarity 0.96 vs new_c++, 0.32 vs c++ |
| MultiStealthyPotential.{h,cpp} | new_c++/CCOptimization/ | only in new_c++ (new main needs it) |
| RepulsiveCCPotential.h | new_c++/CCOptimization/ | |
| DebugSimple.cpp | c++/CCOptimization/ | scratch driver; **excluded** from build |

## apps/CHI/
Source: **`c++/CHI/`** ✅ (no new_c++ copy).

| File | From | Notes |
|------|------|-------|
| main.cpp, ComputationPacking.h | c++/CHI/ | |
| LinealPathFunction.{h,cpp} | c++/CHI/ | `.cpp` byte-identical |
| SpectralDensity.{h,cpp} | c++/CHI/ | |
| VolumeFractionVariance.{h,cpp} | c++/CHI/ | byte-identical |

## apps/CPmanager/  (was legacy `EXC_formats`)
Source: **`new_c++/CPmanager/`** ✅.

| File | From | Notes |
|------|------|-------|
| main.cpp | new_c++/CPmanager/ | ConfigPack ↔ txt manager (similarity 0.97) |

## apps/FFH/
Source: **`c++/FFH/`** ✅ (no new_c++ copy).

| File | From | Notes |
|------|------|-------|
| FFH.{h,cpp} | c++/FFH/ | `.cpp` byte-identical; needs `pairstat` Mtilde + FFTW |
| Tensors.{h,cpp} | c++/FFH/ | `.cpp` byte-identical |
| TwoPhaseDispersion.{h,cpp} | c++/FFH/ | `.cpp` byte-identical |
| main.cpp | c++/FFH/ | |
| FFH_padded.{h,cpp}, FFH_test.cpp, sample1.cpp, sample2.cpp, Dump.hpp | c++/FFH/ | alternates / scratch; **excluded** from build |

## apps/FirstPassageTime/
Source: **`new_c++/FirstPassageTime/`** ✅ (tied with `c++/`; new_c++ authoritative).

| File | From | Notes |
|------|------|-------|
| FirstPassageTime.{h,cpp} | new_c++/FirstPassageTime/ | `.cpp` byte-identical (both trees); needs `nonspherical` |
| main.cpp | new_c++/FirstPassageTime/ | runtime needs its dispersion data files |

## apps/NR/
Source: **`new_c++/NR/`** ✅ (tied with `c++/`; new_c++ authoritative).

| File | From | Notes |
|------|------|-------|
| NR_distribution.{h,cpp} | new_c++/NR/ | |
| examples.cpp | new_c++/NR/ | adapted to lean pairstat API |
| main.cpp | new_c++/NR/ | |
| distance-metrics.cpp | new_c++/NR/ | scratch; **excluded** from build |

## apps/inverse/  (CC-inversion)
Source: **`c++/`** (flat) for the compute sources; entry point **extracted**.
✅ **complete** — builds a runnable `CCInversion.out`, runs (no MathGL/OpenGL), and compiles with **0 warnings in its own sources**; has `README.md` + `run_example.sh`.

| File | From | Notes |
|------|------|-------|
| CCInversion.h | c++/ (flat) | ⚠️ 2026-06-29: holds only the `FrankInvertCCPotential*` family now — `SimpleHarmonic`/`SimpleHarmonicF`/`HS_pot` moved to `misc/misc.h`. |
| CCInversionProblem.{h,cpp} | c++/ (flat) | includes `<LatticeUtil.h>` from misc (2026-06-29) |
| main.cpp | 🆕 | the "CCInversion" command extracted from `c++/Interfaces.h` `cli()`; defines `Verbosity` + `masterDir` |

(`LatticeUtil.h` was promoted out of this app to `misc/` on 2026-06-29 — see the misc section.)

## apps/PairStatCLI/
Source: **`c++/PairStat/`** (driver) + **`new_c++/PairStat/`** (compute).

| File | From | Notes |
|------|------|-------|
| cli.cpp | c++/PairStat/main.cpp | the PairStat CLI driver, renamed (similarity 0.92) |
| Computation.{h,cpp} | new_c++/PairStat/ | newer compute core (similarity 0.96 cpp); this is the `Computation` dropped from the `pairstat` module |

## apps/Statistics/
Source: **`new_c++/Statistics/`** ✅ (no `c++/` copy).

| File | From | Notes |
|------|------|-------|
| gn.{h,cpp} | new_c++/Statistics/ | kernel rewritten (same estimator); dead `__g3_old`/`__g4`/`__g4_old`/`__gn` removed |
| tau_n.{h,cpp} | new_c++/Statistics/ | warning/doc fixes only |
| main.cpp | new_c++/Statistics/ | `gn_1D` flag fix, std::string file names, duplicate core prototypes dropped |
| tests/gen.cpp | new_c++/Statistics/gen.cpp | includes normalized; test-input generator, built by `tests/makefile` |
| tests/check_url_gn.py | new | exact URL g_n validation |
| gn.out (legacy binary), obj/ | — | not copied |

## apps/Tessellation-based/  (tiling / Voronoi)
Source: **`c++/`** (flat); `TilingProblem.cpp` split into three units. ✅ all
objects compile (2026-06-23); has a `README.md`. A runnable binary needs CGAL/voro++/MathGL.

| File | From | Notes |
|------|------|-------|
| Voronoi.{h,cpp} | c++/Voronoi.{h,cpp} | relocated as-is (includes re-normalized) |
| SpacefillingPacking.{h,cpp} | c++/TilingProblem.cpp | 🆕 split: `SpacefillingPacking` class + drivers (legacy lines 16-340) |
| Computation.{h,cpp} | c++/TilingProblem.cpp | 🆕 split: `ImLatGen` + `Computation` hierarchy (legacy lines 345-1362, header-inline); `.cpp` is a stub TU. ⚠️ `ImLatGen` rewritten 2026-06-23 to the current `LatticeGen` API; `MC_nv`→`MC_NV_Ensemble`. |
| main.cpp | c++/TilingProblem.cpp + c++/Interfaces.h | 🆕 `TilingProblemCLI` (legacy 1364-2100) + extracted "Tiling" command. ⚠️ + a local `PoissonPt` generator (absent from modular GenerateConfigs) and `#include <RSA_Gen.h>`/`<HSF_Gen.h>` (2026-06-23). |
| PolyPacking.{h,cpp} | c++/etc1.{cpp,h} | 🆕 extracted `MakeMono2Poly*` family (only tiling consumes it) |

## apps/weighted/
Source: **`new_c++/weighted/`** ✅ (only location).

| File | From | Notes |
|------|------|-------|
| main.cpp, Computation.{h,cpp}, LocalVariance.{h,cpp} | new_c++/weighted/ | byte-identical |
| WeightedStructureFactor.{h,cpp} | new_c++/weighted/ | byte-identical |
| PairCorrelation.{h,cpp} | new_c++/weighted/ | weighted-specific variant (distinct from `pairstat/PairCorrelation`) |

---

## Notes

- **`etc1.{h,cpp}` retired**: no `etc1.*` exists in `modules/`. Its used symbols
  were relocated to `misc/LatticeUtil.h` (`GetLattice_config` /
  `Find_LeastCommonDivider`; first in `apps/inverse/`, promoted to misc 2026-06-29)
  and `apps/Tessellation-based/PolyPacking.{h,cpp}`; the rest were dead/scratch.
- **Makefiles / `config.mk`** are 🆕 (authored for the relocatable modular build),
  not copied from legacy.
- For per-consumer build/blocked status and the migration rules (R1–R7), see
  `PROGRESS.md`; for the tier/layering overview see `README.md`.
