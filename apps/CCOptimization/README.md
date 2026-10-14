# CCOptimization

Ground states and MD samples of stealthy (collective-coordinate) potentials with a soft-core contact
repulsion, in 1–3 dimensions. One binary, `soft_core_stealthy.out`, with two drivers that read their
parameters from standard input:

| Command | Driver | Potential |
|---|---|---|
| `soft_core_stealthy.out [hours [seed [first [verbosity [vtilde]]]]]` | `GetCCO` | One component: `RepulsiveCCPotential` (contact + `ShiftedCCPotential`), direct sums. |
| `soft_core_stealthy.out multi [hours [seed [first [verbosity [vtilde]]]]]` | `GetMultiCCO` | Several species, each with its own constraints: `MultiStealthyPotential`; optional FINUFFT backend. |

`hours` is the time limit (default 144), `vtilde` the shape of v(k): `flat` (default), `overlap` or
`power-law`.

## Build

```sh
make -C apps/CCOptimization                  # Direct sums only
make -C apps/CCOptimization USE_FINUFFT=1    # adds the FINUFFT backend (objects in obj_finufft/)
```

FINUFFT 2.5.1 is found under `$HOME/install`, `/usr/local` or `/usr` (set `FINUFFT_PREFIX`
otherwise). The repository does not compile with `-DNDEBUG`. Bitwise results are tied to the build
flags; the benchmarks below used `-O3 -g` without `-march=native`.

## Input

`GetMultiCCO`, two species of five particles in 2D (the CLI test's input, annotated):

```text
random              initial packing: random | input
2                   dimension
1                   number of different radii
10 1                particle count and relative radius of each radius class
0.01                target packing fraction
0                   relative strength of the contact potential
2                   number of species (collective potentials)
5 0 wavenumber 0 0.1    per species: particle count, S(0), region type (chi | wavenumber), lower, upper
5 0 wavenumber 0 0.1
1                   threads
out/run             save prefix
1                   target number of configurations
ground              mode: ground | MD
finufft             optional keywords, one per line (see below)
run
```

`GetCCO` reads dimension, K1a, K2a, S0, val, sigma, phi, threads, particle number, number of
configurations, `random` or `input` (with a ConfigPack name), the save prefix, the mode, the block
keywords and `run`.

Keywords of the `ground` block (`GetMultiCCO`; `GetCCO` accepts the first three):

| Keyword | Effect |
|---|---|
| `tolerance E` | Energy below which a ground state is accepted (default 1e-14). |
| `maxsteps N` | Evaluation limit of the optimizer (default 10000). |
| `algorithm NAME` | `LBFGS` (default), `LocalGradientDescent`, `ConjugateGradient`, `SteepestDescent`, `MINOP`. |
| `finufft` | Collective densities and, where a cost gate predicts a win, forces with FINUFFT (needs `USE_FINUFFT=1`). |
| `history N` | L-BFGS correction pairs (default 100; 40 was the best measured setting). |
| `combined` | Write the configurations of all species together. |
| `export txt\|binary\|pack\|off` | Output format of accepted configurations. |
| `pairstat Kmax dK Rmax dR samples seed` | Pair statistics of each accepted configuration. |

The `MD` block of `GetMultiCCO` accepts `timestep`, `samplesteps`, `temperature`,
`numequilibration`, `autotimestep`, `restore` and `finufft`.

## Performance notes

- **Direct is the default.** The FINUFFT backend agrees with direct sums to about 1e-14 relative in
  energy and forces. For tiny systems (N of order 100) direct sums stay faster.
- **What the FINUFFT path does.** Densities are Type 1 transforms on persistent plans; forces are
  Type 2 transforms where a deterministic cost gate predicts a win; each species' transforms run as
  particle blocks spread over the threads, and threads without transform work take the contact sum.
  Design, measurements and open items: [docs/finufft/perf/REPORT.md](../../docs/finufft/perf/REPORT.md).
- **Contact cells.** Both drivers size the contact cell list to max(cutoff, (V/N)^(1/d)). Random
  packings start with 10 cells per side, which made the contact sum O(N²) (80% of a 2D N = 10⁴
  FINUFFT evaluation before the fix).
- **Repeatability.** A run with the same build, settings, thread count and input repeats bitwise.
  Endpoints are not bitwise identical to builds before October 2026 (valid ground states either way),
  and direct and FINUFFT endpoints differ because L-BFGS amplifies rounding-level differences.
- **Machine-specific constants.** The cost-gate and block constants were fitted on the development
  machine; refit them on a production node with [tests/finufft_gate/](tests/finufft_gate/).

## Tests and benchmarks

```sh
make -C apps/CCOptimization -j8 all && python3 apps/CCOptimization/tests/test_cli.py
make -C apps/CCOptimization USE_FINUFFT=1 test-finufft test-finufft-forces test-finufft-blocks \
    test-species-concurrent test-indexed-view test-fused-contact
bash potential/tests/finufft_type1/run_tests.sh
make -C apps/CCOptimization USE_FINUFFT=1 bench_finufft.out
./apps/CCOptimization/bench_finufft.out eval|opt DIM N SPECIES CHI THREADS SEED [direct|finufft [HISTORY]]
```
