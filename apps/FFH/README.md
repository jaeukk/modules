# FFH — fast-Fourier-transform homogenization

A C++ implementation of **FFT-based homogenization** for two-phase composite
media, computing the **effective dielectric constant** (permittivity) and the
**effective elastic moduli** (bulk/shear moduli and the full stiffness tensor).
The solver follows the FFT iterative scheme of Moulinec–Suquet with the
accelerated (weak-/strong-contrast) variants of Eyre–Milton.

## Method references

1. D. J. Eyre & G. W. Milton, *A fast numerical scheme for computing the response
   of composites using grid refinement*, **Eur. Phys. J. Appl. Phys. 6**, 41–47 (1999).
2. H. Moulinec & P. Suquet, *A numerical method for computing the overall response
   of nonlinear composites with complex microstructure*, **Comput. Methods Appl.
   Mech. Engrg. 157**, 69–94 (1998).

## Application

Some of the results produced with this code were published in

> J. Kim & S. Torquato, *Multifunctional Composites for Elastic and Electromagnetic
> Wave Propagation*, **Proc. Natl. Acad. Sci. U.S.A. 117**, 8764–8774 (2020).

## How it works

The two-phase medium is discretized onto a regular voxel grid. The Lippmann–Schwinger
equation is solved iteratively: at each stage the polarization field is updated in
direct space, transported to Fourier space via FFTW, convolved with the Green's
function (propagator), and transformed back, until the response field converges.
Both the **weak-contrast** and the **strong-contrast** expansions are available
(`SetMethod` / `SetStrongContrast`).

Problem families (in `FFH.{h,cpp}`):

| class | computes |
|-------|----------|
| `StaticDielectricProblem` | effective permittivity ε_e of a static (electrostatic) two-phase medium |
| `DynamicDielectricProblem` | ε_e including the incident wavevector (electromagnetic wave propagation) |
| `Homo_ElasticProblem_static` | effective bulk/shear moduli and the full stiffness tensor (elastostatics) |
| `DynamicElasticProblem` | elastodynamic effective moduli with a longitudinal/transverse incident wave |

## Layout

| file | contents |
|------|----------|
| [FFH.h](FFH.h) / [.cpp](FFH.cpp) | the solver: `FFH_base` + the dielectric and elastic problem classes, Green's functions, iteration/convergence. |
| [Tensors.h](Tensors.h) / [.cpp](Tensors.cpp) | `ConstantTensor` (a rank-r tensor of complex scalars) and `TensorField` (an FFTW-buffer-backed field over the grid) with inner products and I/O. |
| [TwoPhaseDispersion.h](TwoPhaseDispersion.h) / [.cpp](TwoPhaseDispersion.cpp) | the voxelized two-phase geometry — built from a `SpherePacking` or read from a `.DisArray` file; supports subpixel-smoothing weights. |
| [main.cpp](main.cpp) | the `GetEffEps.out` driver (effective permittivity over an ensemble of `.DisArray` files). Defines `FFH_verbosity`. |

## Dependencies

- **Modules:** `core`, `pairstat` (`Mtilde`).
- **External:** **FFTW3** (with threads), OpenMP, Boost (header-only), GSL (via core).

FFTW is not in the default prefix, so its location is passed at build time.

## Build

```sh
make FFTW_PREFIX=$HOME/miniconda3/envs/pmp        # -> GetEffEps.out
make clean
```

## Driver: `GetEffEps.out`

Computes the effective (static) permittivity of each configuration in an ensemble
and writes the per-configuration values plus the ensemble mean ± standard error.

```sh
./GetEffEps.out <DisArray prefix (printf format)> <num files> <epsilon2> <output prefix> <num threads>
```

| arg | meaning |
|-----|---------|
| `DisArray prefix` | printf-style prefix expanded with the file index to locate each `.DisArray` geometry (e.g. `config_%04d`) |
| `num files` | number of configurations in the ensemble |
| `epsilon2` | permittivity of the inclusion phase (the matrix phase is fixed at ε₁ = 1) |
| `output prefix` | output file prefix |
| `num threads` | OpenMP / FFTW threads |

Each configuration is solved as a `StaticDielectricProblem` (strong-contrast) under a
unit applied field. **Output:** `<output prefix>` — columns `#ID  Re[eps_eff]  Im[eps_eff]`;
if `num files > 1`, the ensemble average and standard error are appended to
`<output prefix>.txt`.

**Example** (100 geometries `dis_0000.DisArray`, …, at inclusion permittivity 10, 8 threads):

```sh
./GetEffEps.out dis_%04d 100 10 results/eps 8
# -> results/eps            (per-config Re/Im effective permittivity)
# -> results/eps.txt        (+ ensemble mean ± SE)
```

The elastic and dynamic problems are exercised through the library API (`FFH.h`)
rather than this driver.

## Scratch / alternate files (not in the default build)

`FFH_padded.{h,cpp}` (zero-padded-grid variant), `FFH_test.cpp` (test `main`),
`sample1.cpp` / `sample2.cpp` (alternate `main`s), and `Dump.hpp` (snippet header)
are kept for reference but are **not** compiled by the makefile.
