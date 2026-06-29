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

Migrated consumers, each with its own makefile. Highlights (see `PROGRESS.md` for
the full per-app status):

- **inverse** — collective-coordinate inversion; **complete**: `make` → runnable
  `CCInversion.out` (needs **no** MathGL/OpenGL), runs, and builds warning-clean.
  Has its own README + example.
- **Tessellation-based** — Voronoi tilings & space-filling packings; all objects
  compile (a *runnable* binary needs CGAL/voro++/MathGL). Has its own README.
- **CHI, NR, FFH, CCOptimization, FirstPassageTime, CPmanager, PairStatCLI** —
  build & run (FFH needs FFTW). `PairStatCLI` is the former pairstat CLI driver.
- **weighted** — weighted S(k) / local variance consumer.

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
