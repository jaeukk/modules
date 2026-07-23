# potential

Pair potentials, dynamics, and energy minimization. Library module (no `main`).

## Contents

| File | Purpose |
|------|---------|
| `Potential.{h,cpp}` | Abstract pair-potential base (energy / force / derivatives). |
| `LJPotential.{h,cpp}` | Lennard-Jones potential. |
| `CollectiveCoordinatePotential.h` | Fixed-cell collective-coordinate (Fourier-space) potentials + the `_ForMC` adapter. |
| `CCPotential_varBox.h` | Variable-cell (box-deformation) CC machinery; the only CC header needing GSL linalg. |
| `CCPotential_Kernels.h` | Catalog of concrete V(k)/S0(k) kernels and the named CC-potential typedefs. |
| `Move.h` | Monte-Carlo move operators (displacement, box scaling). |
| `MonteCarlo.h` | Generic MC drivers (Metropolis, Wang-Landau, parallel tempering). |
| `CoolingSchedule.h` | Annealing schedules (`ExpCool`, `ZeroCool`, adaptive `ThermoCool`). |
| `MC_System.{h,cpp}` | Simulated-annealing MC driver. |
| `MD_System.{h,cpp}` | Molecular dynamics (velocity Verlet, Nosé–Hoover). |
| `StructureOptimization.{h,cpp}` | Energy minimization (NLopt + GSL linear algebra). |
| `AnalyzeData.h` | Autocorrelation / error / correlation-time helpers. |

Header layering (2026-07-07 reorganization; see `CLAUDE.md` for the audit):
`CoolingSchedule.h` → `AnalyzeData.h`; `MonteCarlo.h` → `Move.h` + `CoolingSchedule.h`;
`Move.h` → `Potential.h` only (the MC stack no longer includes the CC potentials);
`CCPotential_Kernels.h` → `CCPotential_varBox.h` → `CollectiveCoordinatePotential.h` → `Potential.h`.

## Dependencies

- **Lower modules:** `core`.
- **External:** GSL, NLopt, OpenMP.

## Build

```sh
make            # compiles all objects (the against-core compile check)
make clean
```

NLopt/GSL prefixes are auto-detected by `../config.mk`.
