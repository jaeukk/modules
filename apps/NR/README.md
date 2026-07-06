# NR — local number distribution $P[N(R)]$

Measures the **local number** $N(R)$ — the number of particle centers inside a
spherical observation window of radius $R$ — over an ensemble of configurations,
and analyzes the resulting distribution $P[N(R)]$ together with its moments
(mean, variance, skewness, excess kurtosis).

For every window radius $R$, **all** sampled $N(R)$ values are stored in a binary
file `*.bin`. Those files are the primary output: they hold the full empirical
histogram of $N(R)$, so the distribution and its higher-order moments can be
(re-)analyzed later without re-scanning the configurations. Standard errors on
the moments are estimated by bootstrap resampling of the stored counts, and each
$P[N(R)]$ can additionally be compared to a Gaussian reference via an $L_2$ norm
and a Kullback–Leibler distance.

## References

Analyses produced with this code appear in:

> 1. M. A. Klatt, J. Kim, T. E. Gartner, and S. Torquato,
>    *Local number fluctuations in ordered and disordered phases of water across
>    temperatures: Higher-order moments and degrees of tetrahedrality*,
>    **J. Chem. Phys.** 160, 204502 (2024).
> 2. S. Torquato, J. Kim, and M. A. Klatt,
>    *Local Number Fluctuations in Hyperuniform and Nonhyperuniform Systems:
>    Higher-Order Moments and Distribution Functions*,
>    **Phys. Rev. X** 11, 021028 (2021).

Please cite the relevant paper when using this code.

## Layout

| file | contents |
|------|----------|
| [NR_distribution.h](NR_distribution.h) / [.cpp](NR_distribution.cpp) | `LocalNumber` (the empirical $N(R)$ histogram: moments, bootstrap SE, binary `WriteBinary`/load); `BootstrapSetup` (what to compute/save); the driver `MC_nv_distribution_PostSampling`; `GetSamplingPoints`; and the Gaussian-reference distance metrics ($L_2$, Kullback–Leibler). |
| [main.cpp](main.cpp) | CLI driver: reads a `ConfigurationPack`, lays out the sampling windows, runs the sweep over $R$, and writes the `.bin` data + analysis. |
| [examples.cpp](examples.cpp) | reference/self-test routines with known distributions (Bernoulli, discrete normal, binomial, Poisson, cubic lattice, RSA). |
| `distance-metrics.cpp` | standalone scratch snippet — **not** part of the build. |

## Dependencies

- **Modules:** `core`, `pairstat` (`NumberVariance`), `models` (`GenerateConfigs`,
  `RandomSequentialAddition`, `HardSphereFluids`).
- **External:** GSL, OpenMP.

All include/library paths come from `../../config.mk`.

## Build

```sh
make            # -> NR.out  (builds core / pairstat / models objects first)
make clean
```

## Run

```sh
./NR.out [printPDF] [printCDF] [periodic] [dx]
```

Optional positional arguments:

| arg | meaning | default |
|-----|---------|---------|
| `printPDF` | `1` = also write the PDF $P[N(R)]$ per radius | `0` |
| `printCDF` | `1` = also write the CDF | `0` |
| `periodic` | `1` = periodic boundary; `0` = hard box (cubic; needs $L > 2R_\max$) | `1` |
| `dx` | window-center spacing: `dx < 0` ⇒ random window centers; `dx ≥ 0` ⇒ a regular grid of centers | `-0.1` (random) |

It then prompts on **stdin** (in order):

```
Configurations to analyze   # ConfigurationPack prefix (the input ensemble)
start index                 # first configuration index to use
number of configurations    # how many to analyze (auto-clamped to what exists)
Random seed for sampling windows
number of sampling windows   # windows per configuration
Maximum radius               # R_max
change in radius             # dR  (radii swept: dR, 2 dR, ..., R_max)
thread number                # OpenMP threads
name of data storage         # prefix for the .bin output files
prefix for the output files  # prefix for the analysis (moments / PDF / CDF) files
```

### Output

- **Binary data** — one file per radius, named `<storage>__<R>.bin`, holding the
  full $N(R)$ histogram for that $R$ across the whole ensemble × windows. Re-run
  with the same storage prefix to reload and re-analyze without rescanning the
  configurations.
- **Analysis** — `<output_prefix>*` text tables: the moments of $P[N(R)]$ vs $R$
  (mean, variance, and, when enabled in `BootstrapSetup`, skewness and excess
  kurtosis), each with a bootstrap standard error, plus the $L_2$ / KL distances
  to a Gaussian; and, if `printPDF`/`printCDF` are set, the distributions
  themselves.

Which quantities are computed/saved is controlled by the `BootstrapSetup` flags
in `main.cpp` (`ComputeVariance`, `ComputeSkewness`, `ComputeKurtosis`,
`ComputeL2norm`, `ComputeKL`, `SaveData`, `PrintPDF`, `PrintCDF`,
`BootstrapSize`).

## Examples

`examples.cpp` provides drivers that build an ensemble with a known $P[N(R)]$
(binomial/Poisson/cubic-lattice/RSA, or a directly sampled Bernoulli / discrete
normal) and run the same analysis pipeline — useful for validating the moment and
distance-metric estimators against analytic results.
