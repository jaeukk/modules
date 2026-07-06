# Tessellation-based — disordered hyperuniform two-phase packings

C++ implementation of the **tessellation-based procedure** for constructing large
realizations of **disordered hyperuniform** two-phase media (packings /
dispersions) from random tessellations.

Starting from a point pattern or monodisperse packing, the method builds its
Voronoi tessellation and grows a particle in each cell to a prescribed fraction
of the **local cell volume**. The resulting **hyperuniform polydisperse packing**
inherits the long-wavelength suppression of density/volume-fraction fluctuations
of the underlying tessellation — i.e. it is hyperuniform, with a spectral density
$\chi_V(k)\to 0$ as $k\to 0$ — while remaining statistically isotropic and
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
| [Voronoi.{h,cpp}](Voronoi.h) | Voronoi-cell computations (cell volume, vertices, neighbours) via **voro++** (`USE_VORO`) / **CGAL** (`USE_CGAL`). |
| [PolyPacking.{h,cpp}](PolyPacking.h) | `MakeMono2Poly*` — grow spheres to a target fraction of each Voronoi-cell volume (the core construction step; needs `USE_CGAL`). |
| [SpacefillingPacking.{h,cpp}](SpacefillingPacking.h) | multiscale space-filling packing (`SpacefillingPacking`, `getSpacefillingPacking*`). |
| [Computation.{h,cpp}](Computation.h) | the `Computation` analysis-task hierarchy: spectral density $\chi_V(k)$, Voronoi cell-volume / vertex distributions, max packing fraction, local number / volume-fraction variance, centroids, and the multiscale family. |
| [main.cpp](main.cpp) | full CLI (`TilingProblemCLI`) — **generates** initial configs internally (RSA / hard-sphere fluids / imperfect lattices / Poisson). |
| [TessellationAnalysisCLI.cpp](TessellationAnalysisCLI.cpp) | the same CLI **decoupled** from the `models` generators (`RSA_Gen`, `HardSphereFluids`, `HSF_Gen`): initial configs are read from files created by other apps. **This is the default build target.** |

## Dependencies

- **Modules:** `core`, `pairstat` (`StructureFactor` / `SpectralDensity` /
  `NumberVariance` / `VolumeFractionVariance`), `models` (the `RSA` base of
  `SpacefillingPacking`, `GenerateConfigs`), `misc` (`Discrepancy` → `MergeSort`).
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

## Usage (TessellationCLI.out)

```sh
./TessellationCLI.out [masterDir]
```

then on **stdin** (top level): `Tiling` enters the CLI, `Exit` quits. Inside, the
two procedures are `VoronoiBased` and `MultiscalePacking`. Typical flow:

```
Tiling
VoronoiBased
  dim 3
  NumThreads 8
  InitialConfiguration LoadPrefix        # read configs made by another app
    <prefix-with-%d>  .txt  <idx_begin>  <idx_end>
  Analyze
    chivComputation ...                  # spectral density chi_V(k)
    VoronoiCellVolumeDistribution ...
    LocalVolumeFractionVariance ...
  SavePrefix <out>
  Start
Exit
```

Add `Construct` (needs `USE_CGAL`) to grow the loaded packings into hyperuniform
polydisperse dispersions before analyzing; `TargetVolumeFraction <phi>` sets the
target. See `TilingProblemCLI` in the source for the full command set.
