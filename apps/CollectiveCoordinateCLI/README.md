# CollectiveCoordinateCLI

A command-line driver for generating and analysing **stealthy hyperuniform**
(and related collective-coordinate) classical ground states in 1, 2 and 3
dimensions. It builds a collective-coordinate potential from a target set of
constrained wavevectors, then finds ground-state point configurations by
gradient optimisation (NLopt L-BFGS, conjugate gradient, MINOP), molecular
dynamics, MD annealing, or Wang–Landau Monte Carlo.

## Authorship & citations

The command interpreter and the collective-coordinate machinery were **written
by Ge Zhang**. This directory is a faithful migration of his legacy
`c++/CollectiveCoordinateCLI.cpp` into the modular `modules/` tree (thin `main`
wrapper + library `-I`/link wiring; see *Migration notes* below). The code was
developed for and used in:

1. Zhang, G., Stillinger, F. H. & Torquato, S. *The perfect glass paradigm:
   Disordered hyperuniform glasses down to absolute zero.* Sci. Rep. **6**,
   36963 (2016).
2. Zhang, G., Stillinger, F. H. & Torquato, S. *Ground states of stealthy
   hyperuniform potentials. II. Stacked-slider phases.* Phys. Rev. E **92**,
   022120 (2015).
3. Zhang, G., Stillinger, F. H. & Torquato, S. *Ground states of stealthy
   hyperuniform potentials: I. Entropically favored configurations.*
   Phys. Rev. E **92**, 022119 (2015).
4. Torquato, S., Zhang, G. & Stillinger, F. H. *Ensemble theory for stealthy
   hyperuniform disordered ground states.* Phys. Rev. X **5**, 021020 (2015).

## Build

```sh
make            # -> CollectiveCoordinateCLI.out
make objs       # compile this app's objects only (no link)
make clean
```

`make` pulls the required library objects from `core`, `pairstat`, `potential`
and `models` (via `include ../../config.mk`, which auto-detects GSL/NLopt). No
module source is modified. External dependencies: **GSL** and **NLopt**.

## Usage

The interpreter is **line-oriented**: each iteration reads one command token,
lets that command consume its argument token(s), then **discards the rest of the
current line**. So put **one command (with its arguments) per line**. A second
command placed after another on the same line is silently dropped — e.g.
`Dimension 2 Exit` sets the dimension and discards `Exit` (the program then just
waits for more input, or ends at end-of-file with status 2). Commands that
*prompt* (e.g. `GeneratePotential` when χ is still 0, or `MultiRun` asking
`Algorithm=`) read that extra input from the following token/line.

Redirect a script (recommended) or type interactively:

```sh
./CollectiveCoordinateCLI.out < commands.txt
```

### Command language

**Parameters (each takes one argument)**

| Command | Meaning |
|---|---|
| `Dimension <d>` | spatial dimension (1/2/3) |
| `NumParticles <N>` (alias `NbrParticles`) | particle count |
| `BoxSize <L>` | cubic box side length |
| `Chi <chi>` | target χ (fraction of constrained degrees of freedom) |
| `NumThreads <n>` (alias `NbrThreads`) | OpenMP threads for the potential |
| `RandomSeed <s>` | RNG seed |
| `CoolingSchedule <v>` (alias `AnnealingSchedule`) | MD-annealing cooling factor |
| `Prefix <p>` / `MDPrefix <p>` | output-file prefix |
| `SampleNumber <k>` / `MDSampleNumber <k>` | number of samples |
| `MDTimeStep <dt>` · `MDTemperature <T>` · `MDStepPerSample <n>` | MD controls |
| `MDEquilibrateSamples <n>` · `MDTimeLimit <sec>` | MD controls |
| `MDAllowRestore <0/1>` · `MDAutoTimeStep <0/1>` · `QuenchAfterMD <0/1>` | MD flags |

**Configuration setup**

| Command | Meaning |
|---|---|
| `GenerateConfiguration` | random config in a cubic box (needs `Dimension`, `NumParticles`, `BoxSize`) |
| `GenerateDensestReciprocalConfiguration` | random config in the densest-reciprocal box (2D/3D) |
| `GenerateFCCConfiguration` | FCC lattice (3D only) |
| `GenerateLatticeConfiguration` | tile a lattice file up to ~`NumParticles` (prompts filename) |
| `ReadConfiguration <file>` | read a `.pos` configuration |
| `LoadConfigPack` | read one entry from a ConfigurationPack (prompts name & index) |
| `ReadProblemFile <file>` | build config + potential from a problem file |
| `AutoRescale` | rescale the config so the stealthy cutoff K = 1 (before `GeneratePotential`) |

**Potential**

| Command | Meaning |
|---|---|
| `GeneratePotential` | build the Type-0 **stealthy** CC potential for the current χ |
| `GenerateOtherPotential` | build another CC potential (prompts χ and a numeric *type*; see the `GeneratePotential` type catalog comment in the source) |

**Optimise / evolve**

| Command | Meaning |
|---|---|
| `CalculateEnergy` | print the current potential energy |
| `Run` | relax via NLopt (L-BFGS) |
| `RunGSL` | relax via GSL conjugate gradient |
| `RunMINOP` | relax via MINOP |
| `RunNewton` | Newton optimisation — **requires `USE_EIGEN`** (see notes) |
| `MultiRun` | `SampleNumber` independent relaxations → ConfigurationPack (prompts `Algorithm`: LBFGS/LocalGradientDescent/ConjugateGradient/SteepestDescent/MINOP) |
| `MD` | NVT molecular dynamics |
| `Anneal` | MD simulated annealing (prompts T_init, T_cut, time step) |
| `WangLandauMonteCarlo` | Wang–Landau MC entropy estimate (prompts several values) |

**Output / control**

| Command | Meaning |
|---|---|
| `Output <prefix>` | write `<prefix>.pos` (and, if plotting were enabled, a figure) |
| `Exit` | free state and quit |

### Example

`run_example.sh` runs a small 2D stealthy ground state end-to-end (N = 40,
χ = 0.3, seed 12345): it generates a random configuration, builds the stealthy
potential, relaxes it to (near) zero energy, writes `stealthy_2D_N40.pos`, then
does a 3-sample `MultiRun` producing `stealthy_2D_N40{,_Success,_InitConfig}`
ConfigurationPacks. The header of that script documents the exact stdin token
order.

```sh
bash run_example.sh
# ... Energy before=10.4776 ... Energy after≈1.7e-21   (a stealthy ground state)
```

## Migration & build notes

- **Figure *rendering* is disabled, but data output is preserved.** The real
  MathGL/OpenGL renderers (`visualization/Plots.cpp`) are unavailable on this box
  (and `Plots.o` cannot compile until `misc/LatticeSumSolver` is migrated — see
  `visualization/CLAUDE.md`). `viz_stubs.cpp` therefore supplies app-local
  implementations of the three referenced symbols:
  - `Plot` and `PlotDirectionalStructureFactor2D` are **no-op** stand-ins (they
    only ever produced bitmap/EPS figures); each prints `[plotting disabled:
    <name>]` to stderr and returns.
  - `PlotFunction_Grace` is the **real legacy implementation copied verbatim**
    from `c++/Plots.cpp` (it is a self-contained plain-text Grace `.agr` *data*
    writer with no MathGL dependency). The Wang–Landau `WangLandauMonteCarlo`
    histogram/entropy/final-entropy `.agr` files are thus written exactly as in
    the legacy program, not dropped.

  Commands that would draw a figure (e.g. `Output`, `Anneal`) still run and still
  write their configuration/data files.
- **The `Debug` command was dropped.** Its dispatch branch and the entire
  `Solvers.h`/`Degeneracy.h` + `CollectiveCoordinateDebug` tail of the legacy
  file were intentionally not migrated (a scratch/experimental sink with no
  stable behaviour). Every other command is preserved verbatim.
- **`USE_EIGEN`.** `RunNewton` and the Hessian-eigenvalue routines need Eigen;
  the default build compiles the non-Eigen stubs (they print an error and
  return). Build with `-DUSE_EIGEN` and Eigen on the include path to enable them.
- **Includes / globals.** Library headers use angle brackets resolved by the
  makefile's `-I` paths (PROGRESS.md R4). The executable owns `Verbosity` and
  seeds the timing globals `ProgramStart`/`TimeLimit` in `main.cpp` (core
  zero-initializes them; leaving `TimeLimit = 0` would make the optimiser's
  time-limit guard exit immediately).
