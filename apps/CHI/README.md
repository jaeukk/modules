# CHI — two-phase-media / sphere-packing statistics command-line interface

An interactive CLI that computes ensemble-averaged **two-phase-media descriptors**
of sphere packings: the spectral density χ_V(k), the lineal-path function L(z),
and the chord-length probability density p(z). It is the **sphere-packing
counterpart of [`apps/PairStatCLI`](../PairStatCLI/)** — same interactive
command-language design and configuration-loader model, but the configurations
are interpreted as a two-phase medium (the spheres are the *particle* phase, their
complement is the *matrix / void* phase) rather than as a bare point pattern, and
the descriptors are those of random heterogeneous media instead of point-process
pair statistics.

> A simplified, packing-oriented variant of the `PairStatisticsCLI` originally
> written by Ge Zhang.

## Provenance

This tool was used to analyze the data published in

> J. Kim & S. Torquato, *Ultradense sphere packings derived from disordered
> stealthy hyperuniform ground states*, **The Journal of Chemical Physics 163**,
> 024902 (2025).

The ultradense packings derived from disordered stealthy hyperuniform ground
states were characterized through their spectral density χ_V(k) and their
lineal-path / chord-length statistics, all computed with this CLI.

## Layout

| file | contents |
|------|----------|
| [main.cpp](main.cpp) | the interactive interpreter (`PairStatisticsCLI`), the configuration-loader helpers (`ReadTextFiles`, `ReadConfigPack`), and `main`. Defines the `Verbosity` global. |
| [ComputationPacking.h](ComputationPacking.h) | the `ComputationPacking` base class and its derived computations (`chivComputation`, `LinealPath`, `ChordLength`, `HoleProbability`). |

The library routines the computations call now live in the `pairstat` module (see
below) — this app ships only its own driver.

## Dependencies

- **Modules:** `core`, `pairstat`. Links the `pairstat` objects
  `SpectralDensity` (χ_V(k)), `StructureFactor` (the sphere form factor
  m̃(k;R), i.e. `Mtilde`), and `LinealPathFunction` (L(z) and p(z), moved out of
  this app into `pairstat`).
- **External:** GSL, Boost (header-only), OpenMP.

All include/library paths come from `../../config.mk`.

## Build

```sh
make            # -> CHI.out  (builds core + pairstat objects first)
make clean
```

## Available computations

Queue a computation by typing its keyword; the CLI then prompts for that
computation's parameters. Type `Available` to list the active ones at runtime, or
`ToDo` for the not-yet-implemented ones.

| command | computes | parameters prompted | output file(s) |
|---------|----------|---------------------|----------------|
| `chivComputation` | isotropic spectral density χ_V(k) | circular `K_max`, linear `K_max`, `K_precision` | `<prefix>_chiv` |
| `LinealPath` | lineal-path function L(z), both phases | `Rmax`, `dR`, `numSamplingLines`, `rng_seed` | `<prefix>_Lz_matrix`, `<prefix>_Lz_particle` |
| `ChordLength` | chord-length density p(z), both phases | `Rmax`, `dR`, `numSamplingLines`, `rng_seed` | `<prefix>_pz_matrix`, `<prefix>_pz_particle` |
| `HoleProbability` | exclusion / hole probability E_V(r) — **stub, not yet implemented** | `Rmax`, `dR`, `numSamplingWindows`, `phase` | `<prefix>_Ev_matrix\|particle` |

`LocalVarianceComputation` (local volume-fraction variance) is present but
commented out — it can be wired to `pairstat`'s `MC_vfv` when needed.

The **phase** convention is `0 = matrix / void`, `1 = particle / spheres`; the
lineal-path and chord-length computations report both phases in one run.

## Interpreter command language

Read on **stdin**; prompts and messages go to stdout. General flow: queue one or
more computations, set the target packing fraction, choose a configuration source,
then run `Calculation`.

| command | effect |
|---------|--------|
| `OutputPrefix <prefix>` | prefix for all output files |
| `NumThreads <n>` | OpenMP threads used by each computation |
| `NumConfig <n>` | number of configurations (auto-clamped by `ReadConfigPack`) |
| `StartIndex <i>` | first configuration index to read |
| `phi2 <φ>` | target packing fraction of the particle phase; **set this before `GetConfigsFunction`** |
| `GetConfigsFunction <type>` | choose the configuration source (see below) |
| `Calculation` | run every queued computation against the source |
| `Available` / `ToDo` | list active / pending computations |
| `Rescale <x>` | *(deprecated)* legacy uniform rescale factor |
| `Exit` | free the queued computations and quit |

### Configuration sources (`GetConfigsFunction <type>`)

| type | reads | notes |
|------|-------|-------|
| `TextFiles` | printf-style filename `prefix`, `start_idx`, `num_configs` | loads each `SpherePacking` (with radii) from a text file `prefix<i>`; if `phi2 > 0`, rescales the radii to the target packing fraction. Auto-sets `NumConfig`. |
| `ReadConfigPack` | pack `prefix` | reads a `ConfigurationPack` of point patterns and turns each into a **monodisperse** packing by choosing the uniform radius that reaches the target `phi2` (**requires `phi2 > 0`**); starts at `StartIndex`. |

## Example session

Compute the spectral density χ_V(k) and the lineal-path function L(z) over a
configuration pack `myPack` of point patterns, interpreted as packings at packing
fraction φ₂ = 0.7, using 8 threads and writing to `results/pk_*`:

```sh
./CHI.out <<'EOF'
OutputPrefix results/pk
NumThreads 8
phi2 0.7
chivComputation
20
40
0.1
LinealPath
5.0
0.02
100000
17
NumConfig 1000
GetConfigsFunction ReadConfigPack
myPack
Calculation
Exit
EOF
```

This queues χ_V(k) (circular `K_max`=20, linear `K_max`=40, bin width 0.1) and
L(z) (`Rmax`=5, `dR`=0.02, 100000 sampling lines, seed 17), opens the pack as
packings at φ₂ = 0.7, runs both over 1000 configurations, and writes
`results/pk_chiv`, `results/pk_Lz_matrix`, and `results/pk_Lz_particle`.

For packings already stored **with radii** as text files
`packing0`, `packing1`, …, read them directly (no `phi2` needed unless you want to
rescale):

```sh
./CHI.out <<'EOF'
OutputPrefix results/cl
NumThreads 8
ChordLength
5.0
0.02
100000
17
GetConfigsFunction TextFiles
packing%d
0
1000
Calculation
Exit
EOF
```
