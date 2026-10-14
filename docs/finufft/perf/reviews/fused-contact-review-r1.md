**BLOCKING — Incremental builds leave incompatible application vtables and crash.**

Locations: `potential/Potential.h:101`, `potential/makefile:34`, `potential/StructureOptimization.cpp:460`.

Trigger: build CollectiveCoordinateCLI at baseline 4586254, apply this change, then run ordinary `make -C apps/CollectiveCoordinateCLI`.

The new dependency rule rebuilds only potential-module objects. Application objects lack header dependencies (`apps/CollectiveCoordinateCLI/makefile:55`; likewise `apps/inverse/makefile:63`). The rebuilt optimizer therefore calls `EnergyAndAllForce` through an old derived-class vtable. Appending the virtual method preserves existing base slots but does not make old derived vtables compatible.

**Evidence run:** I reconstructed the baseline by reversing the supplied diff in a separate disposable copy, built the CLI, applied the diff, and rebuilt without cleaning. I ran:

```sh
./apps/CollectiveCoordinateCLI/CollectiveCoordinateCLI.out <<'EOF'
Dimension 2
NumParticles 40
BoxSize 10
Chi 0.3
NumThreads 4
RandomSeed 12345
GenerateConfiguration
GeneratePotential
CalculateEnergy
Run
CalculateEnergy
Exit
EOF
```

Results:

```text
Baseline:           exit 0; Energy after=1.676225102e-21
Incremental update: return code -11 (SIGSEGV)
Clean rebuild:      exit 0; Energy after=1.676225102e-21
```

The incremental build log confirms that `StructureOptimization.o` rebuilt while `CollectiveCoordinateCLI.o` did not. Ensure affected application objects rebuild when `Potential.h` changes. A successful link alone misses this regression.

**Other validation completed**

| Command | Result |
|---|---|
| `make -C apps/CCOptimization test-fused-contact` | 914 checks passed under both default and passive waiting; 60 Type 2 evaluations each. |
| `make -C apps/CCOptimization test-finufft test-finufft-forces test-species-concurrent test-indexed-view` | All passed, including 168 plan-cache checks and 2,337 concurrency checks per configured run. |
| `bash potential/tests/finufft_type1/run_tests.sh` | All three variants passed: 40,977 / 40,999 / 40,999 checks. |
| `make -C apps/CCOptimization -j8 all && python3 apps/CCOptimization/tests/test_cli.py` | Build and CLI tests passed. |
| `make -C apps/CollectiveCoordinateCLI` and `make -C apps/inverse` | Both linked successfully; incremental runtime failure described above. |

The fused suite covered bitwise agreement, cold AllForce, energy-only behavior, repeated calls, and optimization endpoints.

I additionally built the baseline benchmark and compared:

```sh
bench_finufft.out opt 2 400 3 0.4 4 12345 direct
bench_finufft.out opt 2 400 3 0.4 4 12345 finufft
```

Four interleaved runs per version/backend retained identical final energies, hashes, and counts:

| Backend | Hash | Calls | Baseline median | Changed median |
|---|---|---:|---:|---:|
| Direct | `9427beed288b579f` | 829 | 0.8084 s | 0.7539 s |
| FINUFFT | `d1fa7e5339b16bc6` | 867 | 0.3707 s | 0.3094 s |

The in-process benchmark, `bench_finufft.out eval 2 1000 5 0.3 4 12345`, reported bitwise equality and median evaluation improvements of **7.1% Direct** and **39.9% FINUFFT**. These short measurements support a speedup on the sampled cases; they do not establish the entire claimed performance range.

No additional blocking findings or follow-up findings.