# PairStatCLI — pair-statistics command-line interface

An interactive CLI that computes ensemble-averaged pair / two-point statistics of
point configurations (g(r), S(k), nearest-neighbor and coordination statistics,
local number variance). It is the driver split out of the `pairstat` **library**
module: it ships the high-level `Computation` class hierarchy plus the `cli.cpp`
interpreter, and links against `core` + `pairstat`.

> A simplified version of the `PairStatisticsCLI` originally written by Ge Zhang.

## Layout

| file | contents |
|------|----------|
| [Computation.h](Computation.h) / [.cpp](Computation.cpp) | the `Computation` base class and its derived computations (see below). |
| [cli.cpp](cli.cpp) | the interactive interpreter (`PairStatisticsCLI`) + the configuration-loader helpers, and `main`. Defines the `Verbosity` global. |

## Dependencies

- **Modules:** `core`, `pairstat` (links `PairCorrelation`, `StructureFactor`,
  `NumberVariance`, `ContactStatistics` objects).
- **External:** GSL, Boost (header-only), OpenMP.

All include/library paths come from `../../config.mk`.

## Build

```sh
make            # -> PairStatCLI.out  (builds core + pairstat objects first)
make clean
```

## Available computations

Active (queued by typing the keyword, which then prompts for its parameters):

| command | computes | parameters prompted | output file |
|---------|----------|---------------------|-------------|
| `g2Computation` | isotropic g(r) | `R_max` | `<prefix>_gofr` |
| `Directionalg2Computation` | g(r, θ) | `R_max`, `NumThetaBins` | `<prefix>_x/_y/_g2` |
| `SkComputation` | structure factor S(k) | circular `K_max`, linear `K_max`, `K_precision` | `<prefix>_Sofk` |
| `NearestNeighborStatistics` | nearest-neighbor distance distribution | `SampleDistanceSize` | `<prefix>` files |
| `CoordinationNumber` | coordination-number distribution | `rmin`, `rmax`, `bins` | `<prefix>_Z` |
| `LocalNumberVariance` | local number variance σ²(R) | (reads its own block) | `<prefix>` files |

Type `Available` to list these at runtime, or `ToDo` to list the
not-yet-ported computations (Voronoi-, percolation-, ψ6-based, etc., which are
present but commented out in `cli.cpp`).

## Interpreter command language

Read on **stdin**; output and prompts go to stdout. General flow: queue one or
more computations, choose a configuration source, then run `Calculation`.

| command | effect |
|---------|--------|
| `OutputPrefix <prefix>` | prefix for all output files |
| `NumThreads <n>` | OpenMP threads used by each computation |
| `NumConfig <n>` | number of configurations (auto-set by the ConfigPack loaders) |
| `StartIndex <i>` | first configuration index to read |
| `RandomSampling <0\|1>` | random vs sequential sampling (ReadConfigPack only) |
| `GetConfigsFunction <type>` | choose the configuration source (see below) |
| `Calculation` | run every queued computation against the source |
| `RescaleToUnitDensity` | wrap the source so each configuration is rescaled to unit number density |
| `PrintConfigs` | dump the source configurations as readable text files |
| `AlsoWriteGrace <0\|1>`, `GraceTitle <s>` | also emit Grace plots |
| `Available` / `ToDo` | list active / pending computations |
| `Exit` | free the queued computations and quit |

### Configuration sources (`GetConfigsFunction <type>`)

| type | reads | notes |
|------|-------|-------|
| `PrefixNumberSuffix` | `prefix`, `suffix`, file type | files named `prefix<i>suffix` |
| `PrefixNumberSuffix_withInitialNumber` | + `InitialNumber` | index offset by InitialNumber |
| `ReadConfigPack` | `prefix` | a `ConfigurationPack`; auto-sets `NumConfig` |
| `ReadConfigPack_withInitialNumber` | `prefix`, `InitialNumber` | |
| `SingleConfiguration` | `prefix` | one `.pos` configuration |

File types accepted by the `PrefixNumberSuffix*` loaders: `Pos`, `Xml` (HOOMD),
`Stealth`, `Coordinate`.

## Example session

Compute g(r) up to R=20 and S(k) over a configuration pack named `myConfigs`,
using 4 threads, writing to `results/sample_*`:

```sh
./PairStatCLI.out <<'EOF'
OutputPrefix results/sample
NumThreads 4
g2Computation
20
SkComputation
10
40
0.1
GetConfigsFunction ReadConfigPack
myConfigs
Calculation
Exit
EOF
```

This queues the two computations (each line after the keyword supplies its
prompted parameters), opens the pack (which sets `NumConfig` automatically), runs
both, and writes `results/sample_gofr` and `results/sample_Sofk`.
