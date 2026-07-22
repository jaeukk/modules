# Models CLI

Command-line generator for ensembles of selected disordered point-pattern /
sphere-packing models. Migrated from the legacy `c++/Models/` consumer
(`CLI.cpp`) onto the modular `models` library.

## Build

```sh
make            # builds ./models.out (pulls core + models objects via ../../config.mk)
make clean
```

Depends on `core` and `models` only (GSL + OpenMP).

## Usage

```sh
./models.out -h        # list available models
./models.out models    # interactive generation dialog (reads from stdin)
```

`models` runs a stdin-driven dialog: it reads the model name, that model's
parameters, the number of configurations, an optional input pack (transformation
models only), and the output prefix / type. Output is either a `ConfigPack`
(`<prefix>.ConfigPack` + `.ConfigPackIndex`) or one `SpherePacking` text file per
configuration (`<prefix>__%05d.txt`).

## Available models & input order

Parameters are read in the listed order, one token per line (whitespace works too).

- **Poisson** — spatially uncorrelated (Binomial / Poisson) point pattern.
  ```
  Poisson
  <dimension> <numPrts> <numberDensity> <mode(0=canonical,1=grand-canonical)>
  ```
  (No seed token: the generator seeds from a high-resolution clock.)

- **RSA** — (saturated) random sequential addition of hard hyperspheres.
  ```
  RSA
  <dimension> <numPrts> <numberDensity> <packingFraction> <seed> <numThreads>
  ```
  A saturated packing is generated when `packingFraction >= 0.99 * PhiC(d)` (the
  known saturation value; e.g. pass `1.0`). Otherwise generation stops at the
  requested `packingFraction`. `numberDensity < 0` skips the final rescaling.

- **HSF** — equilibrium hard-sphere fluid sampled by Monte Carlo.
  ```
  HSF
  <dimension> <numPrts> <numberDensity> <packingFraction> <seed> <numThreads> <samplingInterval>
  ```

- **Perturbations** — stochastically displace the particles of a source ensemble
  (requires an input `ConfigPack`).
  ```
  Perturbations
  <dimension> <seed> <distribution> <distribution-parameters...>
  ```
  `distribution` is one of `uniform` (`<maxRadius>`), `gaussian` (`<stddev>`),
  `powerlaw`/`power-law` (`<lengthScale> <exponent>`), or `etc`.

After the model block the dialog reads:

```
<numberOfConfigurations>
<inputPrefix>            # ONLY for Perturbations (source ConfigPack; NumConfig is
                         #   auto-clamped to the pack size)
<outputPrefix>
<outputType>            # ConfigPack | SpherePacking
```

### Note on differences from the legacy CLI

The stdin token order is the **legacy protocol** (all legacy tokens are read in
the legacy positions), mapped onto the reviewed `models` module: `RSA_Gen` (was
`RSA_gen`) has no seed setter — the RSA `<seed>` token is consumed and a warning
notes that RSA seeds from a high-resolution clock; `<numThreads>` is applied via
`SetNumThreads`. HSF applies seed / threads / interval via setters. `HSF_Gen`
uses the argument order `(dimension, numPrts, packingFraction, numberDensity)`
to match `RSA_Gen`.

## Example

`run_example.sh` generates a small RSA `ConfigPack` and a Poisson `SpherePacking`
into the current directory.

## References

- Random sequential addition (RSA):
  G. Zhang and S. Torquato, *Precise algorithm to generate random sequential
  addition of hard hyperspheres at saturation*, **Phys. Rev. E 88, 053312 (2013)**.
