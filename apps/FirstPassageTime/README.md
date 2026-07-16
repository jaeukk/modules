# FirstPassageTime — effective conductivity / diffusion by first-passage-time simulation

C++ implementation of **first-passage-time (FPT)** Brownian-motion ("walk-on-spheres")
simulations that estimate the **effective conductivity** $\sigma_e$ — equivalently
the **effective diffusion coefficient** $D_e$ or dielectric constant — of two-phase
composites made of **disk (2D) / sphere (3D) packings**.

A random walker diffuses through a host phase of conductivity $\sigma_1$ containing
inclusions of conductivity $\sigma_2$. Instead of taking fixed small steps, each
step jumps to the surface of the largest empty concentric sphere ("first-passage"
region), which is exact and lets the walker cover the host efficiently; interface
crossings are handled with the appropriate jump probability set by the conductivity
ratio $\alpha = \sigma_2/\sigma_1$. Averaging the mean square displacement / survival
statistics over many walkers and configurations yields $\sigma_e$.

## References

**Method** (please cite when using this code):

> 1. I. C. Kim and S. Torquato, *Effective conductivity of suspensions of hard
>    spheres by Brownian motion simulation*, **J. Appl. Phys.** 69, 2280–2289 (1991).
> 2. I. C. Kim and S. Torquato, *First-passage-time calculation of the conductivity
>    of continuum models of multiphase composites*, **Phys. Rev. A** 43, 3198–3201
>    (1991). https://doi.org/10.1103/PhysRevA.43.3198
> 3. S. Torquato, I. C. Kim, and D. Cule, *Effective conductivity, dielectric
>    constant, and diffusion coefficient of digitized composite media via
>    first-passage-time equations*, **J. Appl. Phys.** 85, 1560–1571 (1999).
>    https://doi.org/10.1063/1.369287

Used to analyze data in:

> 4. J. Kim and T. M. Truskett, *Geometric model of crack-templated networks for
>    transparent conductive films*, **Appl. Phys. Lett.** 120, 211108 (2022).
>    https://doi.org/10.1063/5.0092946

## Layout

| file | contents |
|------|----------|
| [FirstPassageTime.h](FirstPassageTime.h) / [.cpp](FirstPassageTime.cpp) | `Trajectory` (one walker's state: position, accumulated survival time, displacement) and `FirstPassageTimeParameters` (the walk-on-spheres kernel: single-step jumps, per-configuration and ensemble drivers, 2D and 3D). |
| [main.cpp](main.cpp) | CLI dispatching the run modes below. |

## Dependencies

- **Modules:** `core`, `nonspherical` (`ParticleShapes` / `Dispersions`).
- **External:** GSL, OpenMP.

## Build

```sh
make            # -> FirstPassageTime.out  (builds core + nonspherical objects first)
make clean
```

## Usage

```sh
./FirstPassageTime.out <mode> [args...]
```

`<mode>` = `argv[1]`:

| mode | purpose | arguments |
|------|---------|-----------|
| `2` | **main simulation** — FPT on an ensemble of dispersions from a `ConfigurationPack` | `<ConfigPack prefix> <start idx> <num configs> <save prefix> [log prefix]`, then reads on **stdin**: host σ₁, particle σ₂, target time, #test particles per config, #threads |
| `0` | `Checkerboard` test geometry | `<a> <b> <c> <d>` (floats) |
| `1` | `HoneycombNetwork` test geometry (the crack-templated-network model of ref. 4) | `<5 floats>` |
| `3` | single-trajectory debug on `./config` | `<x> <y>` starting position |
| `4` | validation on a lattice packing | `square|triangular <phi2> <σ₁> <σ₂> …` (reads `./square` / `./triangular`) |

### Typical run (mode 2)

```sh
./FirstPassageTime.out 2 data/spheres 0 100 results/fpt <<'EOF'
1.0        # sigma in host phase
0.0        # sigma in particle phase (0 = perfectly insulating inclusions)
1e4        # target simulation time
10000      # number of test particles per configuration
8          # number of threads
EOF
```

This walks `num-test-particles` walkers through each of the 100 configurations of
the pack `data/spheres`, in parallel, and writes the survival/displacement
statistics (from which $\sigma_e$ / $D_e$ is obtained) under `results/fpt`. An
optional 6th argument is a log prefix for per-trajectory records.

### Other modes (validation / test geometries)

```sh
# mode 0 — checkerboard (args: 4 floats)
./FirstPassageTime.out 0 1.0 0.0 0.5 0.5

# mode 1 — honeycomb network (the crack-templated-network model of ref. 4; 5 floats)
./FirstPassageTime.out 1 0.2 1.0 0.0 0.5 0.5

# mode 4 — periodic lattice validation (reads ./square or ./triangular)
./FirstPassageTime.out 4 square 0.3 1.0 0.0
./FirstPassageTime.out 4 triangular 0.3 1.0 0.0

# mode 3 — single-trajectory debug from position (x, y), reads ./config
./FirstPassageTime.out 3 0.25 0.25
```

> Modes 0/1/3/4 are validation drivers with hard-coded output paths inside
> `main.cpp` (e.g. `./checkerboard-…`); mode 2 is the production driver that reads
> an ensemble and takes its output prefix on the command line. Running with no
> arguments prints a usage line and exits.

### End-to-end example script

```sh
#!/usr/bin/env bash
# run_fpt.sh — effective conductivity of a sphere-packing ensemble
set -euo pipefail
make                                        # build FirstPassageTime.out
mkdir -p results
./FirstPassageTime.out 2 data/spheres 0 100 results/fpt <<'EOF'
1.0
0.0
1e4
10000
8
EOF
echo "done -> results/fpt*"
```
