**No blocking findings. No follow-up findings.** Reviewed only `opt1-r2-fix.diff` and regressions caused by it.

The separate interface resolves the round-1 ABI failure. `Potential.h` matches the reconstructed baseline byte for byte. The fallback in `potential/StructureOptimization.cpp:460–506` preserves the baseline evaluation sequence, and `MultiStealthyPotential` still uses the fused path.

**Incremental-build verification**

I reconstructed the baseline by reversing the supplied diffs, built all three apps, copied in the changed sources, and rebuilt without cleaning:

```sh
make -C apps/CollectiveCoordinateCLI -j8
make -C apps/inverse -j8
make -C apps/CCOptimization -j8
```

Object timestamps and hashes confirmed that **all CollectiveCoordinateCLI and inverse application objects remained unchanged**. CCOptimization rebuilt `main.o`, `MultiStealthyPotential.o`, `MultiComponentPairStat.o`, and `mspb_pack_main.o`.

The supplied CLI input completed successfully:

```text
Baseline:    exit 0; Energy after=1.676225102e-21
Incremental: exit 0; Energy after=1.676225102e-21
```

The inverse executable also exited successfully. Its existing uninitialized time limit causes early termination, so I additionally compiled a test entry point that initializes `TimeLimit`, linking the **unchanged baseline `CCInversionProblem.o`** against the rebuilt optimizer. A six-particle, four-thread BFGS+MINOP reconstruction completed both real and imaginary cases:

```text
run 0
6  6  real  ...  5.96273e-30  ...  BFGS+MINOP
6  6  imag  ...  5.96273e-30  ...  BFGS+MINOP
```

**Commands and results**

| Command | Result |
|---|---|
| `make -C apps/CCOptimization test-fused-contact` | 916 checks passed under both default and passive waiting; 60 Type 2 evaluations each. |
| `make -C apps/CCOptimization test-finufft test-finufft-forces test-species-concurrent test-indexed-view` | All passed: energy oracle, 168 plan-cache checks, all three force runs, 2,337 concurrency checks per run, and indexed-view checks. |
| `python3 apps/CCOptimization/tests/test_cli.py` after incremental rebuild | Passed all CLI checks. |
| `python3 …/check_incremental.py` | Verified unchanged legacy application objects, required CCOptimization rebuilds, and successful application runs. |
| `python3 …/run_sequence.py` | Baseline and fixed objective callbacks produced identical call traces for every switch, with and without gradients. |
| `ownership.out` | Seven checks passed; all six tracked derived destructors ran. |
| `python3 …/run_inverse_entry.py` | Stale inverse driver completed actual optimization; exit 0. |

The callback traces were:

| Switch | Without gradient | With gradient |
|---|---|---|
| 0 | SetConfiguration → Energy | SetConfiguration → Energy → AllForce |
| 1 | SetConfiguration → Energy | SetConfiguration → Energy → basis derivatives |
| 2 | SetConfiguration → Energy | SetConfiguration → Energy → AllForce → basis derivatives |

The ownership harness checked both backends: move construction and assignment preserved bitwise energy/force results, cross-casts worked after moving, and deletion through `Potential*`, `CombinedEnergyForce*`, and `unique_ptr<MultiStealthyPotential>` invoked the derived destructor. Compile-time checks confirmed that copies remain disabled and `Potential` remains abstract, preventing slicing into a base value.