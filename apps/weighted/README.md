# weighted — pair statistics of *weighted* point patterns

An interactive CLI that characterizes the pair statistics — pair correlation
function and structure factor — and other quantities (local variances of the
weights) of **weighted point patterns**: point configurations in which every
particle carries a (complex) weight.

It was written to analyze the data published in

> Torquato, S., Kim, J., Klatt, M. A., Car, R. & Steinhardt, P. J.
> *Hyperuniformity of Weighted Particle Systems.*
> **Phys. Rev. X 16, 011042 (2026).**

A marked point pattern is represented internally as
`MarkedPointPattern = PeriodicCellList<std::complex<double>>`, i.e. an ordinary
point configuration whose per-particle *characteristic* is its weight.

## Usage

**The command language is identical to [`apps/PairStatCLI`](../PairStatCLI/README.md)** —
read commands on stdin, queue one or more computations, choose a configuration
source with `GetConfigsFunction`, then run `Calculation`. See that README for the
shared commands (`OutputPrefix`, `NumThreads`, `NumConfig`, `StartIndex`,
`Calculation`, `Exit`, …).

The **one substantive difference** is the input source: here the configuration
source also needs the per-configuration **weights**, so `ReadConfigPack` takes
three tokens instead of one (see *Input data* below). A few weighted-only knobs
also exist: `ShuffleWeights <0|1>` (randomly permute weights among particles) and
`filter5`/`filter6`/`filter7` (threshold a weight to 0/1; type *before*
`GetConfigsFunction`).

### Build

```sh
make            # -> weighted.out   (builds the core objects first)
make clean
```

Depends on the `core` module only (plus header-only Boost and OpenMP); all
include/library paths come from `../../config.mk`.

### Available computations

Type `Available` at runtime to list these. Each keyword queues a computation and
then reads its parameters from the following tokens:

| command | computes | parameters (in order) | output files |
|---------|----------|-----------------------|--------------|
| `g2Computation` | weighted g₂(r) = ⟨w̄ᵢ wⱼ⟩(r) | `R_max` `dr` | `<prefix>_s2` (or `<prefix>_single-sided-s2`) — rows `(r, Re g₂, Im g₂, err)` |
| `SkComputation` | weighted structure factor S(k) = ⟨\|Σ wⱼ e^{ik·r}\|²⟩/V | `circular_Kmax` `linear_Kmax` `K_precision` | `<prefix>_chif` — rows `(\|k\|, S, err.\|k\|, err.S)` |
| `LocalVariance` | local variance of the weights in a window | `Rmax` `dR` `num_centers` `seed` | `<prefix>_varf[_self]`, `<prefix>_means[_self].txt` |

Per-computation extra options (typed after the computation is queued): `SingleSided`,
`ResolutionPreference`, `ManualBins` (g2); `BraggPeaks`, `CircularKMin`,
`SampleProbability`, `AverageOver…` (Sk); `UseParticleCenters`, `Rmin` (LocalVariance).

## Input data

A marked point pattern is supplied as a **`ConfigurationPack` + companion weights
text files**, wired through:

```
GetConfigsFunction ReadConfigPack  <configpack-prefix>  <weights-prefix>  <column>
```

| token | meaning |
|-------|---------|
| `<configpack-prefix>` | prefix of a `ConfigurationPack`. On disk this is the pair `<prefix>.ConfigPack` + `<prefix>.ConfigPackIndex` (the same binary pack format read by `PairStatCLI`'s `ReadConfigPack`). It also sets `NumConfig` automatically. |
| `<weights-prefix>` | a **printf-style** format string with one integer field for the configuration index — e.g. `…_weights__%04d.txt`. For configuration *i* the loader reads `…_weights__0000.txt`, `…_weights__0001.txt`, … |
| `<column>` | 1-based column of the weights file read as the weight. |

**Weights-file format** — one row per particle, in the **same order** as the
particles in that configuration of the pack:

```
# any line starting with '#' is a comment / column header — skipped
# w1   w2
0.0123   5.0
0.0107   5.0
0.0098   6.0
 …            (NumParticle rows of real, whitespace-separated values)
```

Each weight is parsed as `std::complex<double>`, so a bare real number works (a
complex value may be given in C++ `(re,im)` form). With `<column> = 1` the loader
uses `w1`; with `<column> = 2`, `w2`.

### Specific examples

The two bundled scripts — [`test.sh`](test.sh) and [`test2.sh`](test2.sh) — are
runnable templates (the absolute paths are the author's data directories,
illustrative only). Each pipes a here-string into `weighted.out`; `$config`
expands to the `GetConfigsFunction ReadConfigPack <configpack> <weights-prefix>
<column>` source.

**`test.sh` — weighted structure factor S(k) of a 3D MRJ packing.** Voronoi-cell
volumes as weights (column 1), 50 configurations, 5 threads → `<prefix>_chif.txt`:

```sh
centers="/…/MRJ/3D_MRJ_5000_unitRho"                       # pack: <centers>.ConfigPack(+Index)
weights="/…/MRJ/3D_MRJ_5000_unitRho_weights__%04d.txt 1"   # one weight file per config, column 1
config="GetConfigsFunction ReadConfigPack $centers $weights"

./weighted.out <<< "NumThreads 5
NumConfig 50
$config
SkComputation 3 50 0.1
OutputPrefix $centers
Calculation
Exit
"
```

**`test2.sh` — all three computations on a 2D RSA pattern** (5 configs, weights in
column 1). The runs reproduce the bundled `test1_*`/`test2_*`/`test3_*` output files:

```sh
centers="/…/RSA/test2"
weights="/…/RSA/2D_RSA_0.10_10000_weights__%04d.txt 1"
config="GetConfigsFunction ReadConfigPack $centers $weights"

# test1 — local variance of the weights; 5000 random window centers, seed 123
#         -> test1_varf.txt, test1_means.txt
./weighted.out <<< "NumThreads 2
NumConfig 5
$config
LocalVariance 10 0.05 5000 123
OutputPrefix test1
Calculation
Exit
"
```

The other two runs are the same skeleton with one line changed:

- **test2** — add a `UseParticleCenters` line right after `LocalVariance 10 0.05 5000 123`
  (window centers = particle positions) and use `OutputPrefix test2`
  → `test2_varf_self.txt`, `test2_means_self.txt`.
- **test3** — replace the computation with `g2Computation 10.0 0.05` and use
  `OutputPrefix test3` → `test3_s2.txt`; add a `SingleSided` line for the one-sided
  variant → `test3_single-sided-s2.txt`.

## Layout

| file | contents |
|------|----------|
| [main.cpp](main.cpp) | the interpreter (`PairStatisticsCLI`), the `ReadConfigPackWeights` configuration+weights loader, and `main`. Defines the `Verbosity` global. |
| [Computation.h](Computation.h) / [.cpp](Computation.cpp) | the `Computation` base class and the `g2Computation` / `SkComputation` / `LocalVariance` derived computations. |
| [WeightedStructureFactor.h](WeightedStructureFactor.h) / [.cpp](WeightedStructureFactor.cpp) | weighted S(k) of marked configurations. |
| [PairCorrelation.h](PairCorrelation.h) / [.cpp](PairCorrelation.cpp) | weighted (two-sided & single-sided) g₂(r). |
| [LocalVariance.h](LocalVariance.h) / [.cpp](LocalVariance.cpp) | local variance of the weights over spherical windows. |
| [test.sh](test.sh), [test2.sh](test2.sh) | runnable example scripts (above); [param.txt](param.txt) is a minimal stand-alone LocalVariance parameter file. |

> These `Computation`/`PairCorrelation`/`LocalVariance`/`WeightedStructureFactor`
> files are **distinct** from the `pairstat` module's same-named files: they
> operate on `PeriodicCellList<std::complex<double>>` (weighted) rather than plain
> configurations. Source = the legacy `new_c++/weighted/`.
