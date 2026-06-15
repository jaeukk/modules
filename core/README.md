# core

Foundation module — the types every other module depends on. No internal
dependencies; only external GSL (linear algebra in `GeometryVector`) and Boost
(special functions pulled in via `etc.h`).

## Contents

| File | Purpose |
|------|---------|
| `GeometryVector.{h,cpp}` | Point/vector in up to `MaxDimension` Euclidean dimensions; dot products, GSL-backed matrix helpers. |
| `PeriodicCellList.{h,cpp}` | `PeriodicCellList<T>` — a periodic point pattern with a cell list for O(1) neighbour queries. `Configuration` = `PeriodicCellList<AtomName>`. |
| `RandomGenerator.{h,cpp}` | Reproducible Mersenne-Twister RNG (Boost) with uniform / normal / Poisson draws. |
| `etc.{h,cpp}` | Global constants (`pi`, `MaxDimension`, tolerances) and small shared utilities. |

## Dependencies

- **Lower modules:** none (this is the bottom of the stack).
- **External:** GSL (`-lgsl -lgslcblas`), Boost headers, OpenMP (optional).

## Build

```sh
make            # builds object files + smoke.out
make objs       # object files only (what consumers link against)
./smoke.out     # minimal usage demo (see smoke.cpp)
```

Library paths are resolved by `../config.mk` (auto-detects the GSL prefix). To
override: `make GSL_PREFIX=/path/to/gsl`.

## Minimal usage

See [`smoke.cpp`](smoke.cpp): build a cubic periodic box, insert random points
with a seeded RNG, and count neighbours of a particle with
`IterateThroughNeighbors`.
