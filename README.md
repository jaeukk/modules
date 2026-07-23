# modules/ — modular build of the C++ physics library

A clean, documented, **module-by-module** layout extracted non-destructively from
the legacy flat codebase in `../c++/` (with `../new_c++/` as the newer
authoritative source). The legacy trees are left untouched; consumers migrate
into `modules/` incrementally.

For the running migration status and rules see **`PROGRESS.md`**; for per-file
provenance see **`MODULE_SOURCES.md`**; for build/toolchain notes and conventions
see **`CLAUDE.md`**.

## Dependency tiers

```
          ┌──────────────────────────────────────────────────────────────┐
  apps/   │  inverse  Tessellation-based  CHI  NR  FFH  CCOptimization     │
          │  FirstPassageTime  CPmanager  PairStatCLI  weighted            │
          ├──────────────────────────────────────────────────────────────┤
  leaf    │  phonons   lp (opt Gurobi)   misc   visualization (opt GL/MathGL)
          ├──────────────────────────────────────────────────────────────┤
  mid     │  pairstat      potential      models      nonspherical         │
          ├──────────────────────────────────────────────────────────────┤
  base    │  core   ← everything depends on this                          │
          └──────────────────────────────────────────────────────────────┘
```

Layering rule (R6): a module may only `-I` modules in **lower or equal** tiers.
Library modules live in `modules/<name>/`; migrated application drivers live in
`modules/apps/<name>/`.

## Library modules

| Module | Tier | Purpose | Deps (modules) | External |
|--------|------|---------|----------------|----------|
| [core](core/) | base | GeometryVector, PeriodicCellList, **SpherePacking**, ConfigurationIO, ParameterInput, RandomGenerator, etc; owns `GetKs` | — | GSL, Boost |
| [pairstat](pairstat/) | mid | g(r), S(k), spectral density χ_V, number/volume-fraction variance, contact stats | core | GSL, Boost |
| [potential](potential/) | mid | potentials, MD/MC, structure optimization | core | GSL, NLopt |
| [models](models/) | mid | RSA, hard spheres, random organization, config generators (`RSA_Gen`/`HSF_Gen`) | core (opt pairstat) | GSL |
| [nonspherical](nonspherical/) | mid | ParticleShapes, FormFactors, Dispersions | core | GSL |
| [phonons](phonons/) | leaf | PhononFrequency (+ plotting) | core, potential | GSL, (MathGL) |
| [lp](lp/) | leaf | TorquatoJiao LP jamming (`USE_TJ`) | core | (glpk, Gurobi) |
| [misc](misc/) | leaf | Discrepancy, ParameterSet, Interfaces, LatticeUtil, misc.h (SimpleHarmonic/HS_pot potentials) — header support | core, potential (for `misc.h`) | — |
| [visualization](visualization/) | leaf | Plots, DisplaySpheres, math3d | core | (OpenGL, MathGL) |

Heavy optional libraries (CGAL / voro++ / OpenGL / MathGL / Gurobi) are **off by
default**; the guarded branches compile out so objects build without them.

## Application drivers (`apps/`)

Migrated consumers, each with its own makefile and (unless noted) README.
See `PROGRESS.md` for the full per-app status.

| App | Purpose |
|-----|---------|
| [CCOptimization](apps/CCOptimization/) | Collective-coordinate optimization of stealthy / soft-core hyperuniform ground states (`soft_core_stealthy.out`), including multi-species stealthy potentials with per-species constraints. |
| [CHI](apps/CHI/) | Interactive CLI for ensemble-averaged **two-phase-media descriptors** of sphere packings: spectral density χ_V(k), lineal-path function, and related statistics. |
| [CPmanager](apps/CPmanager/) | Utility to **manage ConfigPack files** (core's `ConfigurationPack` binary ensembles) — convert to/from per-configuration text files. |
| [CollectiveCoordinateCLI](apps/CollectiveCoordinateCLI/) | **Ge Zhang's** driver for generating and analyzing **stealthy hyperuniform** collective-coordinate ground states (MD/MC/annealing + optimizers); used for Sci. Rep. **6**, 36963 (2016), PRE **92**, 022120 (2015), PRE **92**, 022119 (2015), PRX **5**, 021020 (2015). |
| [FFH](apps/FFH/) | **FFT-based homogenization** of two-phase composites: effective dielectric constant / elastic response (needs FFTW). |
| [FirstPassageTime](apps/FirstPassageTime/) | **First-passage-time** (walk-on-spheres) Brownian-motion simulations estimating effective conductivity of dispersions. |
| [Models](apps/Models/) | Command-line generator (`models.out`) for ensembles of disordered models: Poisson, RSA, equilibrium hard spheres, perturbed lattices (legacy stdin protocol). |
| [NR](apps/NR/) | Local number statistics N(R) in spherical observation windows over an ensemble — number-variance analysis. |
| [OrderMetric_PeriodicNetwork](apps/OrderMetric_PeriodicNetwork/) | Two programs for **nonspherical** two-phase media: orientationally-averaged χ_V(k) of dispersions (`GetChivDispersions.out`) and surface-area coefficient / volume-fraction variance of periodic networks (`NetworkBV.out`). |
| [PairStatCLI](apps/PairStatCLI/) | Interactive CLI for ensemble-averaged pair statistics of point configurations: g₂(r), S(k), nearest-neighbor and coordination statistics (the former `pairstat` CLI driver). |
| [Tessellation-based](apps/Tessellation-based/) | **Tessellation-based procedure** for constructing large disordered hyperuniform two-phase media; Voronoi/space-filling analysis CLI (`TessellationCLI.out`, needs `USE_VORO=1`). |
| [inverse](apps/inverse/) | Collective-coordinate **inversion** (`CCInversion.out`) — complete, warning-clean, no MathGL/OpenGL needed. |
| [polygon_packing](apps/polygon_packing/) | Placeholder — currently empty (no sources migrated yet). |
| [weighted](apps/weighted/) | Pair statistics of **weighted/marked** point patterns: weighted S(k), two-/one-sided g₂, local variance of marks (PRX **16**, 011042 (2026)). |

## Build

**Toolchain:** there is no native Windows g++/make here — build through **WSL**
(repo at `/mnt/e/30_Codes/modules`). Every makefile does `include config.mk`,
which **auto-detects** the GSL/NLopt prefix across machines (`~/install`,
`~/lib/gsl`, `/usr/local`, `/usr`). Override with `make GSL_PREFIX=/path`.

```sh
make -C core                    # foundation (+ smoke.out)
make -C pairstat objs           # library objects (the CLI lives in apps/PairStatCLI)
make -C potential objs          # library objects
make -C models objs             # library objects
make -C examples                # end-to-end demo consumer (g2_sk_demo.out)
make -C apps/inverse            # -> runnable CCInversion.out
make -C apps/Tessellation-based # compile objects (binary needs CGAL/voro++/MathGL)
```

Optional features are command-line toggles handled by `config.mk`:
`USE_OPENMP` (default 1), `USE_CGAL`, `USE_VORO`, `USE_TJ`, `USE_MATHGL`, and
per-module flags like `WITH_SK`.

## Provenance & migration

- `MODULE_SOURCES.md` — which legacy copy each file came from, and why.
- `PROGRESS.md` — migration rules (R1–R7) and per-module / per-consumer status.
- `CLAUDE.md` — build/toolchain how-to, conventions, recurring fix patterns,
  dated change log.
- `examples/` — the reference consumer and migration recipe for porting the
  remaining subprojects off their broken `../` / absolute-path makefiles.
