# FINUFFT cost model: calibration

With the `finufft` backend, `MultiStealthyPotential::AllForce` (and `ShiftedCCPotential::AllForce`) computes a
species' collective-coordinate forces with FINUFFT Type 2 instead of the direct per-particle loop
when a cost model predicts a win by `FinufftForceMargin` (1.3); `ShiftedCCPotential::FinufftForceGate` decides
jointly over the species that run concurrently, and also how many particle blocks each Type 2 runs as. Type 1
(rho) is split into particle blocks by `FinufftBlocks`, from the inputs only. The model's constants
(`FinufftCostCall`, `FinufftCostParticle`, `FinufftCostMode`, `FinufftCostPair` for Type 2 and the direct loop,
`FinufftRhoCostCall`, `FinufftRhoCostParticle`, `FinufftRhoCostMode` for Type 1, in
`potential/CollectiveCoordinatePotential.h`) are machine-specific: they only choose the faster path and the work
split, never the accuracy (the Type 1 block count does change rho in the last bits, deterministically). Refit them
on a new CPU, compiler or FINUFFT version, e.g. the production node. Keep all output outside the repository (e.g.
the project's data area).

From `apps/CCOptimization`, on an otherwise idle node:

    make finufft-gate                                   # builds calibrate_gate.out and check_gate.out
    ./calibrate_gate.out > $DATA/gate_calibration.txt   # 142 cases, several minutes; or DIM N CHI [SKEW] for one
    python3 tests/finufft_gate/fit_gate.py $DATA/gate_calibration.txt

`fit_gate.py` prints the seven constant lines (tightest linear upper envelopes of the Type 2 times and of the
Type 1 block times per dimension, smallest long-loop direct pair rate) and a margin check: accepted cases must be
measured wins. Paste the lines over the constants (and update the date and machine in the `ponytail:` comments above
them), then:

    make finufft-gate && ./check_gate.out        # gated vs Type-2-disabled AllForce; ratio < 1 where accepted
    make test-finufft-forces                     # Combined expects designed gate decisions at T=4 and 8
    make test-finufft-blocks                     # expects several Type 1 blocks in its large cases

`FinufftBlockRatio` and `FinufftBlockMinWork` (same header) set how far the transforms split; they were chosen with
a thread-CPU schedule model of whole evaluations (2D/3D, N = 2000 and 10000, S = 2 and 5, T = 4 and 8).

If `test-finufft-forces` reports "gate decision differs from the expected one", the new constants moved a case
across the threshold: update that case's expected decisions in `Combined` (tests/test_finufft_forces.cpp).
