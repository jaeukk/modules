# modules/ — modular build of the C++ physics library

A clean, documented, **module-by-module** layout extracted non-destructively from
the legacy flat codebase in `../`. The legacy top-level files and the existing
application subdirectories are left untouched; consumers migrate to `modules/`
incrementally.

## Dependency tiers

```
            ┌─────────────────────────────────────────────────────────┐
  T-leaf    │  viz (opt: OpenGL/MathGL)                                │
            ├─────────────────────────────────────────────────────────┤
  T-high    │  inverse (opt: Gurobi)      tiling (opt: CGAL/voro++)    │
            ├─────────────────────────────────────────────────────────┤
  T-mid     │  pairstat      potential      models                     │
            ├─────────────────────────────────────────────────────────┤
  T-base    │  core   ← everything depends on this                     │
            └─────────────────────────────────────────────────────────┘
```

Layering rule: a module may only `-I` modules in **lower or equal** tiers.

| Module | Tier | Purpose | Deps (modules) | External |
|--------|------|---------|----------------|----------|
| [core](core/) | base | GeometryVector, PeriodicCellList, RandomGenerator, etc | — | GSL, Boost |
| [pairstat](pairstat/) | mid | g(r), S(k), number variance | core | GSL, Boost |
| [potential](potential/) | mid | potentials, MD/MC, optimization | core | GSL, NLopt |
| [models](models/) | mid | RSA, hard spheres, random organization | core (opt pairstat) | GSL |
| inverse | high | inverse problems, CC inversion, tilings-LP | core, potential, pairstat | GSL, NLopt, (Gurobi) |
| tiling | high | Voronoi, tiling problems | core, pairstat | CGAL, voro++ |
| viz | leaf | DisplaySpheres, Plots | core | OpenGL, MathGL |

`inverse`, `tiling`, `viz` are scaffolded incrementally; the heavy optional
libraries (CGAL / voro++ / OpenGL / MathGL / Gurobi) are off by default.

## Build

Every module/consumer makefile does `include ../config.mk`, which **auto-detects**
the GSL/NLopt prefix across machines (HOME `~/lib/gsl`, WSL `~/install`, DELLA
`/usr/local`). Override with `make GSL_PREFIX=/path NLOPT_PREFIX=/path`.

```sh
make -C core            # foundation + smoke.out
make -C pairstat        # + pairstat.out CLI
make -C potential       # library objects
make -C models          # library objects  (make WITH_SK=1 to enable S(k) features)
make -C examples        # end-to-end demo consumer (g2_sk_demo.out)
```

Optional features are command-line toggles handled by `config.mk`:
`USE_OPENMP` (default 1), `USE_CGAL`, `USE_VORO`, `USE_TJ`.

## Provenance & migration

- `MODULE_SOURCES.md` — which legacy copy each file came from, and why.
- `examples/` — the reference consumer + the step-by-step migration recipe for
  porting the remaining subprojects (CHI, FFH, NR, HSF, …) off their broken
  `../` / absolute-path makefiles onto `modules/`.
