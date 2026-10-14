# FINUFFT speedup of the multi-component optimization

October 14, 2026. This report supersedes the performance part of [../REPORT.md](../REPORT.md)
(September 28), where only the densities rho(k) used FINUFFT and forces stayed direct.

## Result

With the integrated design, one L-BFGS evaluation of `GetMultiCCO` (2D, N=2000, five species of 400,
packing fraction 0.2, 4 OpenMP threads) costs **1.45 ms** at chi_s=0.4, versus 11.65 ms with Direct
sums and 7.21 ms with the September FINUFFT path. End to end, that ground-state optimization takes
**15 s instead of 73 s (Direct) or 44–54 s (old FINUFFT)**. Energies and forces agree with Direct to
1e-14 relative, and every FINUFFT run is bitwise repeatable across processes.

**How to use it.** Build with `make -C apps/CCOptimization USE_FINUFFT=1` (FINUFFT is detected under
`$HOME/install`, `/usr/local` or `/usr`; set `FINUFFT_PREFIX` otherwise), then add the keyword
`finufft` to the `ground` or `MD` parameter block before `run`. Direct stays the default. FINUFFT's
own thread count stays 1; parallelism comes from the `# threads` input. For tiny systems (N of order
100) keep Direct: FINUFFT is no slower than before there, but still slower than Direct.

## Where the time went (2D, N=2000, S=5, 4 threads, ms per evaluation)

| Part | chi_s=0.1 Direct | chi_s=0.4 Direct | chi_s=0.4 old FINUFFT |
|---|---|---|---|
| Configuration build + `SetConfiguration` (S sub-configuration rebuilds) | 0.45 | 0.5 | 0.5 |
| rho(k), five species one after another | 1.05 | 4.2 | 1.15 |
| CC + contact energy | 0.15 | 0.16 | 0.16 |
| `AllForce` (direct O(N M) CC force sum) | 1.26 | 4.6 | 4.6 |

The old FINUFFT rho cost was almost entirely plan creation (238 us per species call against 44 us
for a reused plan), species ran one after another on single-threaded FINUFFT calls, and the force sum
was untouched. The app did not link FINUFFT at all, so production runs used Direct everywhere.

## Process

1. Shared plumbing (c352d19, dd2fb45): `USE_FINUFFT=1` build toggle with a separate object directory
   (ODR) and header dependencies, the `finufft` keyword, `MultiStealthyPotential::SetRhoBackend`, the
   Direct-vs-FINUFFT oracle `tests/test_finufft_multi.cpp`, and the benchmark `tests/bench_finufft.*`.
   It also fixed the mode-snapping guard of `GetRho_FinufftType1`, which rejected every lattice mode
   with index above about 43/d at eps=1e-12 and so fell back to Direct silently in 1D and at large N.
   A Codex review of the plumbing found no regression; its compressed-box force finding was rejected
   with evidence (Direct's own error against the exact zero force is the same size) and turned into
   a scale-aware near-zero bound, and a NaN gap in the oracle was fixed.
2. Brainstorm with Codex ([brainstorm-codex.md](brainstorm-codex.md)): four distinct mechanisms, with
   Type 3 transforms and an automatic Direct/FINUFFT crossover deferred.
3. Each option implemented on its own branch from c352d19 and audited by a Codex critic (GPT-6 Astra),
   up to three rounds, under the agreed rules: a finding blocks only with a realistic trigger and a
   wrong-physics, crash, lost-repeatability or slowdown impact; later rounds review only the fix diff;
   stop after a round with no blocking findings.
4. Integration of all four (user's choice), audited the same way.

## The four options

| Option | Mechanism | Audit rounds | Blocking findings and fix |
|---|---|---|---|
| A, Type 2 forces | CC forces as one Type 2 transform per species (ntransf = d) | 2 | Small systems got slower (3–12x): a deterministic cost gate keeps Type 2 only where it is predicted to win |
| B, persistent plans | One guru Type 1 plan and mode map per species, bitwise fingerprint for invalidation | 1 | None (a test target that skipped the plan-cache test was fixed) |
| C, species concurrency | One parallel region; FINUFFT species run concurrently on one thread each | 2 | A FINUFFT species that falls back to Direct ran single-threaded (5x): fallback species now keep the full team |
| D, indexed views | Species potentials read the full configuration through an index list instead of S rebuilt sub-configurations | 1 | None |

Standalone, FINUFFT backend, median ms per evaluation (ratio to the old FINUFFT path):

| Case | Old FINUFFT | A | B | C | D |
|---|---|---|---|---|---|
| 2D chi_s=0.1, T=4 | 3.27 | 3.20 (0.98) | 2.40 (0.73) | 2.35 (0.72) | 3.08 (0.94) |
| 2D chi_s=0.4, T=4 | 6.45 | 3.45 (0.53) | 5.65 (0.88) | 5.59 (0.87) | 6.20 (0.96) |
| 2D chi_s=0.4, T=8 | 7.28 | 5.23 (0.72) | 3.97 (0.54) | 3.81 (0.52) | 6.50 (0.89) |
| 3D chi_s=0.4, T=4 | 11.55 | 11.64 (1.01) | 10.38 (0.90) | 9.26 (0.80) | 10.83 (0.94) |

A attacks the force sum, B and C the rho/energy part, and D the per-evaluation setup, so they combine.

## Integrated design

One cache per species holds the bitwise input fingerprint, the eligibility result and mode map, a
Type 1 plan and a Type 2 plan (ntransf = d), and the particle angles, which are marked as belonging
to the cached rho. Plans are made only in a serial prepare step before each parallel region.
`Energy` runs FINUFFT species concurrently, one thread each, while free threads start the contact sum.
`AllForce` runs, per FINUFFT species, rho and (if the gate accepts) Type 2 concurrently, then one
particle loop adds the contact force. Every Direct or declined species keeps the c352d19 arithmetic
bitwise. The Type 2 gate is a deterministic joint rule over the concurrent species (never timing),
with constants calibrated on this machine; see `apps/CCOptimization/tests/finufft_gate/` to refit
them on another node. Integration audit: round 1 found wrong Type 2 forces after copying an evaluated
potential and a Direct-path contact-schedule slowdown at small N; both fixed; round 2 found nothing
blocking. Its follow-up (contact chunk too coarse at N=1e4) was applied as a one-line cap.

## Integrated benchmark (interleaved with the base, another heavy job on the machine)

Median ms per evaluation:

| Case | Base Direct | Base FINUFFT | Integrated FINUFFT | vs Direct | vs base FINUFFT |
|---|---|---|---|---|---|
| 2D N=2000 chi_s=0.1 T=4 | 3.59 | 3.58 | 1.44 | 2.5x | 2.5x |
| 2D N=2000 chi_s=0.4 T=4 | 11.65 | 7.21 | 1.45 | 8.0x | 5.0x |
| 2D N=2000 chi_s=0.1 T=8 | 2.08 | 2.89 | 0.76 | 2.8x | 3.8x |
| 2D N=2000 chi_s=0.4 T=8 | 9.98 | 7.48 | 0.77 | 12.9x | 9.7x |
| 3D N=2000 chi_s=0.1 T=4 | 5.38 | 6.72 | 4.02 | 1.3x | 1.7x |
| 3D N=2000 chi_s=0.4 T=4 | 15.18 | 12.17 | 5.07 | 3.0x | 2.4x |
| 2D N=10000 chi_s=0.4 T=4 | 223.1 | 136.2 | 18.1 | 12.3x | 7.5x |
| 3D N=100 chi_s=0.4 T=4 | 0.17 | 1.74 | 0.31 | 0.6x | 5.7x |

Full ground-state optimizations (three-pass L-BFGS, T=4), seconds:

| Case | Base Direct | Base FINUFFT (2 runs) | Integrated FINUFFT (2 runs) | Final energy (Direct) |
|---|---|---|---|---|
| 2D chi_s=0.1 | 9.8 | 7.4 / 6.9 | 3.9 / 3.9 | 8e-20 |
| 2D chi_s=0.4 | 72.7 | 43.7 / 53.8 | 15.0 / 14.7 | 5e-20 |
| 3D chi_s=0.1 | 1.3 | 1.7 / 1.7 | 1.0 / 1.0 | 7e-20 |
| 3D chi_s=0.4 | 19.0 | 16.6 / 15.2 | 7.1 / 8.5 | 3e-19 |

Both integrated runs of each case ended at bitwise-identical coordinates. Evaluation counts differ
slightly between backends (for example 3701 against 3448 calls), as expected from rounding-level
evaluator differences; the September report documents that L-BFGS amplifies these.

## Correctness evidence

`make -C apps/CCOptimization test-finufft test-finufft-forces test-species-concurrent test-indexed-view`
(Direct-vs-FINUFFT oracle in 1–3D and skew cells, persistent-plan cache and invalidation with
injected FINUFFT failures, Type 2 gate decisions and forces, bitwise repeatability at 4 and 8
threads and under passive waiting and nesting, cold `AllForce` in MD order, indexed views against
explicit sub-configurations), `bash potential/tests/finufft_type1/run_tests.sh` (all three builds),
the repeatability API tests and probe of `potential/tests/finufft_repeatability`, an
AddressSanitizer build of the plan-cache test, and the default build with `tests/test_cli.py`.
Direct-backend energies and forces are bitwise identical to dd2fb45 in 1–3D at 1, 4 and 8 threads.

## Caveats found after the benchmarks

**Build flags.** All repository code in these benchmarks was built with the `config.mk` defaults,
`-std=c++17 -g -O3 -fopenmp`: no `-march=native` and no `-DNDEBUG`. Without `NDEBUG`,
`GeometryVector` records its dimension and asserts in every operation, which the direct force loop
executes once per particle-mode pair. The FINUFFT library is a native (AVX-512) build. The speedups
therefore lean in FINUFFT's favor; base and new used the same flags, so the ranking holds, but the
factors would shrink somewhat in a release build. The repository does not currently compile with
`-DNDEBUG`: `core/PeriodicCellList.h` (two `displacement_periodic` overloads),
`pairstat/NumberVariance.h` and `pairstat/VolumeFractionVariance.cpp` read `GeometryVector::Dimension`,
which exists only in debug builds.

**Ground states are not bit-identical to the old build.** Each run is deterministic for a fixed build,
flags, CPU, settings and thread count, so endpoint hashes are a valid comparison. From the same initial
state (seed 1, N=2000, five species, 2D/3D, chi_s 0.1/0.4):

| Variant | Same endpoint as the base FINUFFT path |
|---|---|
| B (persistent plans), D (indexed views) | yes, bitwise, in all four cases |
| A (Type 2 forces) | yes where the gate keeps the direct loop; no at 2D chi_s=0.4 |
| C (species concurrency) | no in three cases (summation order of the species energy changed); yes at 3D chi_s=0.4 |
| Integrated | no in all four cases |

All endpoints are ground states (total energy 2e-20 to 4e-19 by Direct re-evaluation). The differences
are not caused by the stopping tolerance: trajectories diverge early (call counts differ by tens to
hundreds), because a last-bit change in the energy flips near-tie line-search decisions and the
stealthy ground-state manifold is degenerate for chi < 0.5. Rerunning a seed from an existing dataset
with the new build will not reproduce its configuration. Whether old and new produce statistically
equivalent ensembles (g2, S(k) outside the constraints, nearest-neighbor distances over many seeds) is
not yet checked.

A visual explainer of the algorithmic changes is in [explainer.html](explainer.html).

## Open items

1. The Type 2 gate constants are machine-specific (Ryzen 5 7500F, WSL2, FINUFFT 2.5.1). Refit them
   on the Nurion nodes before production use there; the expected decisions in the Combined test of
   `test_finufft_forces.cpp` follow the constants.
2. For tiny systems FINUFFT rho is still slower than Direct; an automatic per-species crossover for
   Type 1 (the deferred brainstorm idea) would remove the need to choose by hand.
3. A failed Type 2 plan or execution is visible only as `GetLastForceBackendUsed() == Direct`.
4. `SetFinufftNumThreads(>1)` is not guarded inside the concurrent region (not reachable from the CLI).
5. `potential/tests/contact_reproducibility/run_finufft_threads.py` still expects the historical
   incoming request of 4 FINUFFT threads and fails on fresh builds (predates this work).
6. DUCC0 starts an idle thread pool on the first FINUFFT call, as before this work.
7. `ShiftedCCPotential_S2_MultiConfig` with FINUFFT selected builds an unused Type 1 plan in `AllForce`.
8. Direct-vs-FINUFFT endpoint identity is still not claimed (unchanged from September).

## Provenance

Branches (local, not pushed): `jaeukk/cco-optA-type2-forces` (6e57558, 093a5bb),
`jaeukk/cco-optB-persistent-plans` (5f91e7d, 79b5371, 9f0bb99), `jaeukk/cco-optC-species-concurrent`
(6ce4b8f, 7917475), `jaeukk/cco-optD-indexed-view` (95a608c), and `jaeukk/cco-finufft-integration`
(merges 68b8abf, 4e93bad, a01f36c, 5c11101; fixes 7f5933a to 5d209ca; gate calibration c8147ff). Codex reviews are in
[reviews/](reviews/); task contracts and dispatch records are in `~/.claude/multiagent/tasks/cco-finufft-*`.
