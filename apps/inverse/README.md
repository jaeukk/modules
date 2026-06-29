# inverse — Collective-coordinate inversion

> **Status: complete.** `make` builds a runnable `CCInversion.out` (no MathGL/OpenGL),
> the example runs, and the app's own sources compile **warning-clean** (`-Wall`).

Reconstructs point configurations from a prescribed subset of their **collective
coordinates** (constraints on the structure factor / Fourier-space density), and
studies how many distinct real-space configurations are consistent with a given
set of constraints — i.e. the (non-)uniqueness of the inverse map.

The reconstruction is posed as energy minimization: a `FrankInvertCCPotential`
objective (in [CCInversion.h](CCInversion.h)) penalizes the mismatch between a
trial configuration's collective coordinates and the target's, and is driven to
its ground states from many random starts by the optimizers in the `potential`
module (BFGS+MINOP, steepest descent, conjugate gradient, local gradient
descent).

## Reference

This work is published in:

> J. Kim, G. Zhang, F. H. Stillinger, and S. Torquato,
> *Inversion problems for Fourier transforms of particle distributions*,
> **J. Stat. Mech.: Theory Exp.** 2018, 113302 (2018).

Please cite it when using this code.

## Layout

| file | contents |
|------|----------|
| [CCInversion.h](CCInversion.h) | CC-inversion potentials (`FrankInvertCCPotential*`). |
| [CCInversionProblem.h](CCInversionProblem.h) / [.cpp](CCInversionProblem.cpp) | uniqueness-test drivers (`FrankInvCCUniqueTest*`) and the CLI entry `CCInversionProblem1`. |
| [main.cpp](main.cpp) | entry point — the `CCInversion` command extracted from the legacy `Interfaces.h` dispatcher. |

> The generic potentials `SimpleHarmonic`/`SimpleHarmonicF`/`HS_pot` and the
> lattice helpers `GetLattice_config`/`Find_LeastCommonDivider` (`LatticeUtil.h`)
> were moved into the **`misc`** module on 2026-06-29; this app now includes
> `<LatticeUtil.h>` from there.

## Dependencies

- **Modules:** `core`, `pairstat` (`StructureFactor` → `GetKs`), `potential`
  (`Potential`, `StructureOptimization`), and `misc` (`Discrepancy` → `MergeSort`,
  compiled into this app's `obj/` by the makefile; `LatticeUtil.h`).
- **External:** GSL, NLopt, OpenMP. No MathGL/OpenGL needed — the `CCInversion`
  command does not touch the visualization implementations.

All include/library paths come from `../../config.mk`.

## Build

```sh
make            # -> CCInversion.out
make clean
```

## Run

```sh
./CCInversion.out [outputBaseDir]      # outputBaseDir defaults to "."
```

It then reads a small command language on **stdin**:

```
CCInversion <printAllSolutions 0|1>
<NKSkip> <dim> <minN> <maxN> <dN> <numConfigs> <Nc> <systType1> <systType2>
<methodNum1> <methodNum2> <EnergyTolerance> <delta1> <delta2> <d_delta>
<M1> <M2> <numThreads> <doProxMetric>
Exit
```

Parameters (see `CCInversionProblem1` in [CCInversionProblem.cpp](CCInversionProblem.cpp)):

| token | meaning |
|-------|---------|
| `NKSkip` | number of central k-points to omit from the constraint set |
| `dim` | spatial dimension |
| `minN`..`maxN` step `dN` | target system sizes (number of particles) |
| `numConfigs` | number of random target configurations per case |
| `Nc` | number of random initial guesses per target |
| `systType1..2` | target ensemble: `0` = weakly perturbed lattice (amplitude `delta`), `1` = Poisson, `2` = non-Bravais (1D only) |
| `methodNum1..2` | optimizer: `0` = BFGS+MINOP, `1` = steepest descent, `2` = conjugate gradient, `3` = local gradient descent |
| `EnergyTolerance` | ground-state energy threshold |
| `delta1..delta2` step `d_delta` | perturbation amplitude sweep (for `systType` 0) |
| `M1`, `M2` | number-of-constraints range, as `N+M1 .. N+M2` |
| `numThreads` | OpenMP threads |
| `doProxMetric` | `1` to also compute solution-proximity metrics |

### Output

`CCInversionProblem1` writes one table per `(N, optimizer, system type)` to a
**CWD-relative** path that it does **not** create automatically:

```
results/Reconst/dE_<EnergyTolerance>/Nk_skip_<NKSkip>/<dim>D<systName>__<method>__N_<N>.txt
```

so create that directory before running (the example script does this for you).
Each table reports, per number of constraints `M`, the spread of distinct
solutions found (`maxSols`/`minSols`/`av.NumSols`), the probability of finding
multiple solutions, the worst ground-state energy, and the largest position
error.

## Example

[run_example.sh](run_example.sh) builds the binary if needed, creates the output
directory, and runs a tiny 1D weakly-perturbed-lattice case (`N = 4`,
`delta = 0.05`, BFGS+MINOP):

```sh
./run_example.sh
```

It prints where the result table was written and shows its first lines. Edit the
parameter block in the script (e.g. raise `maxN`, `numConfigs`, `Nc`, or sweep
`delta`) to reproduce larger studies from the paper.
