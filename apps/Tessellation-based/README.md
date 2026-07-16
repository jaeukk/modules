# Tessellation-based — disordered hyperuniform two-phase packings

C++ implementation of the **tessellation-based procedure** for constructing large
realizations of **disordered hyperuniform** two-phase media (packings /
dispersions) from random tessellations.

Starting from a point pattern or monodisperse packing, the method builds its
Voronoi tessellation and grows a particle in each cell to a prescribed fraction
of the **local cell volume**. The resulting **hyperuniform polydisperse packing**
inherits the long-wavelength suppression of density/volume-fraction fluctuations
of the underlying tessellation — i.e. it is hyperuniform, with a spectral density
$\tilde{\chi}_V(k)\to 0$ as $k\to 0$ — while remaining statistically isotropic and
disordered. A complementary *multiscale* space-filling construction fills space
with an ever-finer sequence of particle sizes.

## References

**Method** (please cite when using this code):

> 1. J. Kim and S. Torquato, *Methodology to construct large realizations of
>    perfectly hyperuniform disordered packings*, **Phys. Rev. E** 99, 052141 (2019).
> 2. J. Kim and S. Torquato, *New tessellation-based procedure to design perfectly
>    hyperuniform disordered dispersions for materials discovery*,
>    **Acta Mater.** 168, 143–151 (2019).

The hyperuniform polydisperse packings produced here were also used in:

> 3. J. Kim and S. Torquato, *Theoretical prediction of the effective dynamic
>    dielectric constant of disordered hyperuniform anisotropic composites beyond
>    the long-wavelength regime*, **Opt. Mater. Express** 14, 194 (2024).
> 4. J. Kim and S. Torquato, *Extraordinary Optical and Transport Properties of
>    Disordered Stealthy Hyperuniform Two-Phase Media*,
>    **J. Phys.: Condens. Matter** 36, 225701 (2024).
> 5. J. Kim and S. Torquato, *Effective electromagnetic wave properties of
>    disordered stealthy hyperuniform layered media beyond the quasistatic
>    regime*, **Optica** 10, 965–972 (2023).
> 6. S. Torquato and J. Kim, *Nonlocal Effective Electromagnetic Wave
>    Characteristics of Composite Media: Beyond the Quasistatic Regime*,
>    **Phys. Rev. X** 11, 021002 (2021).

## Layout

| file | contents |
|------|----------|
| [Voronoi.{h,cpp}](Voronoi.h) | Voronoi-cell computations (cell volume, vertices, neighbours) via **voro++** (`USE_VORO`) / **CGAL** (`USE_CGAL`). Carries a file-local `MergeSort` / `Merge_Sort` helper (copied from the former `misc/Discrepancy`) used to order particles along a coordinate in the 1-D branch — so the app needs no `misc` dependency. |
| [PolyPacking.{h,cpp}](PolyPacking.h) | `MakeMono2Poly*` — grow spheres to a target fraction of each Voronoi-cell volume (the core construction step; needs `USE_CGAL`). |
| [SpacefillingPacking.{h,cpp}](SpacefillingPacking.h) | multiscale space-filling packing (`SpacefillingPacking`, `getSpacefillingPacking*`). |
| [Computation.{h,cpp}](Computation.h) | the `Computation` analysis-task hierarchy: spectral density $\chi_V(k)$, Voronoi cell-volume / vertex distributions, max packing fraction, local number / volume-fraction variance, centroids, and the multiscale family. |
| [main.cpp](main.cpp) | full CLI (`TilingProblemCLI`) — **generates** initial configs internally (RSA / hard-sphere fluids / imperfect lattices / Poisson). |
| [TessellationAnalysisCLI.cpp](TessellationAnalysisCLI.cpp) | the same CLI **decoupled** from the `models` generators (`RSA_Gen`, `HardSphereFluids`, `HSF_Gen`): initial configs are read from files created by other apps. **This is the default build target.** |

## Dependencies

- **Modules:** `core`, `pairstat` (`StructureFactor` / `SpectralDensity` /
  `NumberVariance` / `VolumeFractionVariance`), `models` (the `RSA` base of
  `SpacefillingPacking`, `GenerateConfigs`). *(No `misc` dependency: the one
  merge-sort helper the Voronoi code needs is copied file-locally into
  `Voronoi.cpp`.)*
- **External:** GSL; and the optional heavy libraries below.

### Optional heavy libraries (feature flags)

| flag | library | enables |
|------|---------|---------|
| `USE_VORO` | **voro++** (`voro++-dev`) | Voronoi-cell analyses (the tessellation itself) |
| `USE_CGAL` | **CGAL** | the polydisperse **construction** (`MakeMono2Poly`) |
| `USE_MATHGL` | **MathGL** | on-the-fly plots |

`config.mk` auto-detects a voro++ install (`VORO_PREFIX`, searching
`$HOME/install`, `/usr/local`, `/usr` for `include/voro++/voro++.hh`) and emits
`VORO_INC` / `VORO_LIB`. Without these libraries the corresponding branches are
`#ifdef`-guarded out and the objects still compile.

## CLI drivers

Two interactive drivers share the same `Tiling` command tree (`VoronoiBased` /
`MultiscalePacking`); they differ only in how initial configurations enter:

| driver | binary | initial configurations | built by default? |
|--------|--------|-------------------------|-------------------|
| **[TessellationAnalysisCLI.cpp](TessellationAnalysisCLI.cpp)** — the **analysis-only CLI** | `TessellationCLI.out` | **loaded from files** produced by other apps (`InitialConfiguration LoadPrefix …`) | ✅ **yes — this is the makefile target** |
| [main.cpp](main.cpp) — the generator-coupled CLI | (`Tiling.out`) | **generated internally** (RSA / hard-sphere fluids / imperfect lattices / Poisson) | ❌ no — needs the unmigrated `models/HSF_Gen.h` |

So the binary the makefile produces, **`TessellationCLI.out`, is the analysis-only
CLI**: it never generates configurations itself — you feed it packings / point
patterns created by another app and it runs the Voronoi-based and multiscale
`Computation` analyses (optionally the `USE_CGAL` polydisperse construction
first). This mirrors the file-scanning workflow of `apps/PairStatCLI`.

## Build

```sh
make objs               # compile the library objects (no heavy libs needed)
make USE_VORO=1         # -> TessellationCLI.out : the config-reading CLI, with
                        #    voro++ Voronoi analyses enabled
make USE_VORO=1 USE_CGAL=1   # additionally enables the MakeMono2Poly construction
make clean
```

`main.cpp` (the generator-coupled CLI) additionally needs `models/HSF_Gen.h`,
which has not been migrated into the modular tree yet, so it is **not** part of
the default build — use `TessellationAnalysisCLI.cpp` and feed it configurations
produced by another app.

## Running the analysis-only CLI (`TessellationCLI.out`)

**Important:** the binary takes **no analysis arguments on the command line** — it
reads an interactive **command script from `stdin`**. The optional positional
argument only sets the working directory:

```sh
./TessellationCLI.out [masterDir]        # masterDir (argv[1]) = working dir; default "./"
```

So you drive it by piping a script in (or typing the commands interactively):

```sh
./TessellationCLI.out < script.txt
# or a heredoc:
./TessellationCLI.out <<'EOF'
Tiling
VoronoiBased
  dim 3
  NumThreads 8
  NumConfigs 100
  InitialConfiguration LoadPrefix
    config_%04d  .txt  0  100
  Analyze
  SavePrefix result
  chivComputation                 0.0  10.0  0.1  0
  LocalVolumeFractionVariance     5.0  0.1  2000
  VoronoiCellVolumeDistribution   10.0 0.05 0
  Start
Exit
Exit
EOF
```

What this does: loads the 100 files `config_0000.txt … config_0099.txt`, runs the
spectral density $\chi_V(k)$, the local volume-fraction variance, and the Voronoi
cell-volume histogram, and writes outputs prefixed with `result`.

### Command tree

- Top level: `Tiling` enters the CLI, `Verbosity <n>` sets log level, `Exit` quits
  the program. Inside `Tiling` the two procedures are **`VoronoiBased`** and
  **`MultiscalePacking`**; each is ended by `Start` (runs) then `Exit` (leaves).
- `VoronoiBased` — whitespace-separated tokens, any order, then `Start`:
  | token | arguments | meaning |
  |-------|-----------|---------|
  | `dim` | `d` | spatial dimension |
  | `NumThreads` | `n` | OpenMP threads |
  | `NumConfigs` | `N` | ensemble size |
  | `InitialConfiguration LoadPrefix` | `prefix .type idx_begin idx_end` | load configs: `prefix` is a printf pattern with `%d`, `.type` is `.txt`, indices run over the half-open range `idx_begin` to `idx_end`. **The files must already exist** — if any is missing the loader re-prompts and eats the following tokens. |
  | `Analyze` | — | enable the analysis pass (required for the tasks to run) |
  | `SavePrefix` | `prefix` | output filename prefix |
  | `Construct` | — | grow the loaded packings into hyperuniform polydisperse dispersions first (**needs `USE_CGAL`**); `TargetVolumeFraction <phi>` sets the target |
  | `Start` | — | run the configured tasks over the ensemble |
- Analysis **tasks** (add any subset before `Start`; each reads its own parameters):
  | task | parameters |
  |------|-----------|
  | `chivComputation` | `CircularKmax LinearKmax Kprecision lowStorage(1/0)` — isotropic $\chi_V(k)$ → `<prefix>_chiv` |
  | `LocalVolumeFractionVariance` | `Rmax dR numSamplingPts` → `<prefix>_localVolVar` |
  | `VoronoiCellVolumeDistribution` | `maxCellVolume binWidth analyzeInitPacking(1/0)` |
  | `VoronoiCellVertexDistribution` | (cell-vertex histogram) |
  | `LocalNumberVarianceComputation` | (local number variance) |
  | `CentroidComputation` | `outputEachConfig(1/0)` |

See `TilingProblemCLI()` in [TessellationAnalysisCLI.cpp](TessellationAnalysisCLI.cpp)
for the complete grammar.
