# potential

Pair potentials, dynamics, and energy minimization. Library module (no `main`).

## Contents

| File | Purpose |
|------|---------|
| `Potential.{h,cpp}` | Abstract pair-potential base (energy / force / derivatives). |
| `LJPotential.{h,cpp}` | Lennard-Jones potential. |
| `CollectiveCoordinatePotential.h` | Collective-coordinate (Fourier-space) potential. |
| `Move.h` | Monte-Carlo move operators (displacement, box scaling, CC moves). |
| `MonteCarlo.h` | Generic MC framework with cooling schedules. |
| `MC_System.{h,cpp}` | Simulated-annealing MC driver. |
| `MD_System.{h,cpp}` | Molecular dynamics (velocity Verlet, Nosé–Hoover). |
| `StructureOptimization.{h,cpp}` | Energy minimization (NLopt + GSL linear algebra). |
| `AnalyzeData.h` | Autocorrelation / error / correlation-time helpers. |

## Dependencies

- **Lower modules:** `core`.
- **External:** GSL, NLopt, OpenMP.

## Build

```sh
make            # compiles all objects (the against-core compile check)
make clean
```

NLopt/GSL prefixes are auto-detected by `../config.mk`.
