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

## Follow-up: fused contact pass and L-BFGS history

An independent Codex review of the result proposed two more improvements
([reviews/codex-next-improvements.md](reviews/codex-next-improvements.md)).

**Fused contact energy and force (landed).** A gradient evaluation now walks each particle's
neighbor list once for both its contact-energy slot and its contact force, and the species' Type 2
work overlaps the contact walk instead of running in its own phase. The combined call lives on a
separate interface, `potential/CombinedEnergyForce.h`, which `MultiStealthyPotential` implements and
`Objective` finds with `dynamic_cast`. `Potential`'s vtable is unchanged, so applications whose
objects predate the change keep working after an incremental build (the first audit round caught a
crash from an earlier version that added a virtual to `Potential`). Results are bitwise identical to
calling `Energy()` then `AllForce()`, and full optimizations end at the same coordinates as before.
One evaluation got 26–54% faster with FINUFFT (2D N=2000: 1.13 to 0.64 ms median at chi_s=0.4,
T=4; N=10000: 17.6 to 8.5 ms) and 3–13% faster with Direct. Two Codex rounds: one blocking finding
(the vtable change), fixed; round 2 clean.

**Configurable L-BFGS history (landed, default unchanged).** A time split of whole optimizations
showed L-BFGS's own work at 24–38% of the time in the FINUFFT cases (configuration rebuild and
gradient conversion only 2–3%); with the fused evaluator, L-BFGS at history 100 costs more per call
than the potential (1.09 against 0.61 ms at 2D chi_s=0.4). `StructureOptimization_LBFGSHistory`
(default 100) and the `history N` keyword in the `ground` block make it adjustable. Minimum of three
interleaved full optimizations (seed 1, 4 threads, fused build):

| Case | history 100 | history 40 | history 20 |
|---|---|---|---|
| 2D chi_s=0.1 | 1.52 s | 0.98 s (0.64x) | 0.82 s (0.54x) |
| 2D chi_s=0.4 | 6.59 s | 4.94 s (0.75x) | 5.46 s (0.83x) |
| 3D chi_s=0.4 | 5.55 s | 4.20 s (0.76x) | 4.46 s (0.80x) |

A shorter history needs more evaluations (three seeds, history 40 against 100: +0% to +25%; history
20 at 2D chi_s=0.4: +61%). Every history from 10 to 100 reached valid ground states in all 12 runs
and repeated bitwise. A shorter history changes trajectories and endpoints for both backends, so the
default stays 100; `history 40` is the measured best general setting.

## Follow-up: contact cells and the thread-allocation measurement

The question was whether threads sit idle while FINUFFT species run one per thread, so that the
stealthy part should be split further. An instrumented copy of `EnergyAndAllForce` (2D/3D,
N=2000/10000, S=2/5, T=4/8, chi_s 0.4) found the team only 1–4% idle, but the contact walk at 80% of an
evaluation at 2D N=10000 and growing as N^2. Cause: random packings start from `GetUnitCubicBox(d, 0.1)`,
whose cell list keeps 10 cells per side through `Rescale`, and `PrepareIterateThroughNeighbors` never
re-bins, so each contact query scanned about 0.09 N particles.

**Contact cells (commit 24727d9).** `MultiStealthyPotential::PrepareContactCells` re-bins the
per-evaluation copy to cells of width max(rc, (V/N)^(1/d)) in `SetConfiguration`, `SetSecureRadii` and
`SetPacking`. FINUFFT evaluation, 2D, S=5, T=4: N=2000 1.17 to 0.77 ms; N=10000 18.3 to 3.66 ms. 3D at
N=2000 already had cells near 1.7 rc and is unchanged; 3D N=10000 contact time fell 25–35%. The width
rule matched or beat a fixed width of 0.55. Contact sums now visit neighbors in a different order, so 2D
endpoints change in the last bits; they remain ground states and repeat bitwise across processes (the 3D
N=2000 endpoint hash is unchanged).

**GetCCO contact cells (commit 930d18a).** The single-component driver had the same coarse cells.
`RepulsiveCCPotential` points at the caller's configuration instead of owning a copy, so re-binning
inside it would cost a copy per evaluation (15–20% slower at 2D N=500). Instead `GetCCO` wraps
`GetInitConfigs` to set the cell size to max(sigma, (V/N)^(1/d)) once per configuration; the optimizer's
per-evaluation copies (`GetConfig` copies `origConfig`) and MD's working copy inherit it. With FINUFFT rho
(`RepulsiveCCPotential` through the API; 2D, chi 0.4, T=4, interleaved): one evaluation N=2000 1.6 to
1.2 ms, N=10000 about 17 to 3.9 ms, no change at N <= 1000; a full N=10000 optimization 34 to 20 s. 3D is
within noise up to N=5000 (cells already near the cutoff) and about 4% faster at N=10000. The CLI's
`GetCCO` runs Direct sums only, where the gain is about 3% at 2D N=5000. 2D endpoints change in the last
bits, remain ground states and repeat bitwise; `test_cli.py` runs a GetCCO ground state twice and compares
the packs byte for byte.

**Thread allocation after the fix.** With the contact walk cheap, a thread-CPU schedule model (each task
timed with `CLOCK_THREAD_CPUTIME_ID`, list-scheduled on T threads; robust to the machine's load) bounded
what a perfect split of the species work could gain: 2D S=2 T=8 2.35–2.67x; 3D S=2 T=4 1.42–1.57x;
3D S=2 T=8 about 3x; 3D S=5 T=8 1.42–1.57x; S=5 T=4, and generally S >= T with equal species, about
1.0x. Ordering species largest first gains nothing for equal species. The hypothesis held for S < T
(binary mixtures, 3D) but not for S >= T, which motivated the particle blocks below.

## Follow-up: species work in particle blocks

With one FINUFFT species per thread, a run with fewer species than threads left most of the team idle
while each species' transforms ran (one 3D species of 2000-5000 particles at chi 0.4 takes 4-20 ms on one
thread, mostly spreading and interpolation). Each species' transforms now run as particle blocks.

**Design.** A block is a contiguous range of the species' local particle order (block b of B covers
[bN/B, (b+1)N/B)) with its own persistent guru plan; plans are still made only in the serial prepare step.

- Type 1 (rho): B1 blocks, each transform giving its own mode grid; rho at a constraint is the sum of the
  blocks' modes in block order. B1 = `FinufftBlocks(1)` depends only on N, G, d and eps: the largest B for
  which a block's predicted spreading work stays at least twice the work every block repeats (call, FFT and
  deconvolution on the mode grid) and at least 100 us, at most 16. So rho is bitwise identical for every
  thread count; B1 = 1 (every small species, and every species of the chi 0.4 benchmarks below) is bitwise
  the previous path, and B1 > 1 changes rho in the last bits (within FINUFFT tolerance of the single
  transform).
- Type 2 (forces): blocks of the same kind. Every block deconvolves and transforms the same coefficients and
  FINUFFT interpolates each point on its own, so the forces are bitwise those of the unsplit transform for
  any block count (checked against the simple interface in 1-3D). The count therefore may depend on the
  team: `FinufftForceGate` picks one count b <= floor(T/F) for its F species (each species also capped by
  `FinufftBlocks(2)`) that minimizes the predicted region time, max(longest species chain, total work / T),
  with the contact walk as other work (`MultiStealthyPotential::ContactCostParticle`), and keeps b > 1 only
  if that wins by 1.2 (`FinufftSplitMargin`). Its acceptance rule now costs the blocks (work = blocks x block
  time, wall time added estimated as ceil(blocks/T) x the longest block).
- Schedule: one region per evaluation, one monotonic dynamic loop over the Type 1 blocks (largest predicted
  time first), then the Type 2 blocks of the species with several. The thread that finishes a species' last
  Type 1 block (seq_cst countdown) continues with it: rho reduction in block order (or the direct
  fallback), energy, Type 2 input, and either its single Type 2 block or the release of its several blocks;
  a thread holding a block of a species not yet released takes contact chunks meanwhile. Then the remaining
  contact chunks and one barrier, as before; `Energy`, `AllForce` and `EnergyAndAllForce` share it, so the
  fused call stays bitwise equal to the two calls. An earlier version with a barrier between a Type 1 and a
  Type 2 phase cost 3-8% on small 2D systems under load and was replaced.
- A lone `ShiftedCCPotential` runs its blocks in parallel over `ParallelNumThread` in `GetRho` and
  `TryFinufftAllForce`. Direct and fallback species keep their previous path.

The Type 1 cost constants come from the extended calibration tool (`calibrate_gate.out` now also times one
Type 1 block; `fit_gate.py` fits both envelopes); the Type 2 constants are unchanged.

**Thread-CPU schedule model** (every task timed serially with `CLOCK_THREAD_CPUTIME_ID`, minimum over 7
evaluations, then list-scheduled on T threads as each build schedules it; best of two passes; chi 0.4;
ms per fused evaluation; blocks = Type 1 / Type 2 tasks). Cases with identical work in both builds scatter
by up to 8%, which is the noise floor of this machine under load:

| Case | Base | Blocks | Speedup | Blocks (T1/T2) |
|---|---|---|---|---|
| 2D N=2000 S=2 T=4 | 0.41 | 0.35 | 1.19 | 2/0 |
| 2D N=2000 S=2 T=8 | 0.33 | 0.35 | 0.95 | 2/0 |
| 2D N=2000 S=5 T=4 | 0.23 | 0.25 | 0.96 | 5/0 |
| 2D N=2000 S=5 T=8 | 0.12 | 0.12 | 1.02 | 5/0 |
| 2D N=10000 S=2 T=4 | 1.13 | 1.15 | 0.98 | 2/0 |
| 2D N=10000 S=2 T=8 | 1.09 | 0.83 | 1.31 | 2/8 |
| 2D N=10000 S=5 T=4 | 1.12 | 1.10 | 1.02 | 5/0 |
| 2D N=10000 S=5 T=8 | 0.56 | 0.56 | 1.02 | 5/0 |
| 3D N=2000 S=2 T=4 | 3.01 | 2.20 | 1.37 | 2/4 |
| 3D N=2000 S=2 T=8 | 3.82 | 2.44 | 1.57 | 2/8 |
| 3D N=2000 S=5 T=4 | 3.77 | 3.51 | 1.07 | 5/0 |
| 3D N=2000 S=5 T=8 | 1.90 | 1.91 | 1.00 | 5/0 |
| 3D N=10000 S=2 T=4 | 14.66 | 11.13 | 1.32 | 2/4 |
| 3D N=10000 S=2 T=8 | 17.29 | 8.28 | 2.09 | 2/8 |
| 3D N=10000 S=5 T=4 | 11.67 | 11.64 | 1.00 | 5/0 |
| 3D N=10000 S=5 T=8 | 6.13 | 6.02 | 1.02 | 5/0 |

**Wall clock** (`bench_finufft.out eval ... finufft`, three interleaved runs per build, load 3-7 on the 8
threads; median of the per-run medians, and minimum, of SetConfiguration + EnergyAndAllForce in ms). Medians
at T=8 are dominated by preemption; the minimums are the steadier comparison:

| Case | Base median / min | Blocks median / min | Ratio of minimums |
|---|---|---|---|
| 2D N=2000 S=2 T=4 | 0.47 / 0.34 | 0.47 / 0.34 | 0.98 |
| 2D N=2000 S=5 T=4 | 0.44 / 0.36 | 0.44 / 0.35 | 1.01 |
| 2D N=10000 S=2 T=8 | 2.67 / 1.55 | 3.47 / 1.56 | 1.00 |
| 2D N=10000 S=5 T=8 | 1.46 / 1.16 | 1.49 / 1.18 | 0.98 |
| 3D N=2000 S=2 T=4 | 4.00 / 2.99 | 2.99 / 2.28 | 1.31 |
| 3D N=2000 S=2 T=8 | 9.83 / 3.03 | 7.71 / 2.37 | 1.28 |
| 3D N=10000 S=2 T=4 | 19.52 / 16.02 | 15.39 / 12.51 | 1.28 |
| 3D N=10000 S=2 T=8 | 24.17 / 17.87 | 22.57 / 14.94 | 1.20 |
| 3D N=10000 S=5 T=4 | 14.23 / 12.72 | 14.24 / 12.55 | 1.01 |

All 16 cases of the model table were run; the omitted ones are within 0.98-1.08 by minimum. The split
3D cases gain 1.2-1.3x here against 1.3-2.1x in the model: the other jobs leave fewer free cores than
the model's idle team. The 2D N=10000 S=2 T=8 gain of the model (1.31x) does not show under this load.

Full optimizations (`bench_finufft.out opt`, seed 1, history 100, two interleaved processes per build):

| Case | Base (s) | Blocks (s) | Calls | Final energy (FINUFFT / Direct) | Hash |
|---|---|---|---|---|---|
| 3D N=2000 S=2 T=8 | 3.35 / 3.34 | 2.56 / 2.72 | 324 in all four | 3.589e-20 / 3.597e-20 | 8783c1a687821511 in all four |
| 2D N=10000 S=2 T=8 | 30.3 / 30.3 | 31.7 / 30.6 | 1996 in all four | 1.239e-20 / 1.285e-20 | 7ba17c859dee174b in all four |

Both builds end at bitwise the same coordinates: at chi 0.4 every Type 1 count is 1 and Type 2 blocks change
no bit. The 2D N=10000 optimization is dominated by L-BFGS (15 ms per call against about 2 ms per
evaluation), so the evaluation change is within its noise.

**Checks.** New `test_finufft_blocks` (`make test-finufft-blocks`): rho is bitwise the block-order sum of
simple-interface Type 1 transforms of the blocks and within tolerance of the single transform and of
Direct, bitwise identical for 1, 4 and 8 threads and on repetition; Type 2 forces in several blocks are
bitwise those of the unsplit simple-interface transform; mixed species (split, unsplit, FINUFFT fallback,
Direct) at T = 1, 4, 8 with 2, 4 and 6 Type 2 blocks give bitwise identical results, every Type 1 block runs
once per evaluation, the cold `AllForce` equals the warm one and `EnergyAndAllForce` equals `Energy` then
`AllForce`; block plans are trimmed and completed with the particle count, remade on a tolerance change,
never shared by copies, released by assignment, and injected makeplan and setpts failures in one block fall
back cleanly. Mutations caught: reducing rho in thread-completion order (T-dependence and repetition checks
fail), reversing the block order (block-order check fails), and a one-particle offset in the Type 2 blocks.
Existing suites pass with two adapted observations: `test_species_concurrent` counts `GetRho` cache misses,
which FINUFFT species no longer make (their rho comes from block tasks; the lazy-fill check stays for
`GetRho`); and in `test_finufft_forces` the mid-size 3D case now uses Type 2 with eight threads (its blocks
run in parallel), so the thread-count check moved to a dense case too small to split, and the mid-size case
checks one block at T=1 against several at T=8, bitwise. AddressSanitizer runs of the block and plan-cache
tests are clean.

**Audit.** One Codex round over 854eea9..930d18a (contact cells, particle blocks, GetCCO cells;
[reviews/blocks-review-r1.md](reviews/blocks-review-r1.md)) found no blocking issue. Besides the suites
above it ran the block test with `OMP_THREAD_LIMIT=1`, with `OMP_THREAD_LIMIT=2 OMP_DYNAMIC=true`, and inside
an outer two-thread region (nested levels 1 and 4) without deadlock, and repeated the 3D N=2000 S=2
optimization hash. Its one follow-up, a stale open item on GetCCO cells, is fixed.

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
9. Particle blocks: the Type 1 constants, `FinufftBlockRatio`, `FinufftSplitMargin` and
   `MultiStealthyPotential::ContactCostParticle` are fitted on this machine under load; refit with
   `tests/finufft_gate` (and re-time the contact walk) on the production node. The contact estimate ignores
   the packing fraction.
10. Particle blocks: the model gains of the split cases (up to 2.1x) were only partly visible on this loaded
    machine (1.2-1.3x); confirm on an idle node. A species with B1 > 1 changes rho in the last bits against
    earlier builds (none of the chi 0.4 benchmarks does; low chi and large N do).
11. Particle blocks reach only part of the schedule-model bound: 3D S=5 T=8 stays at 1.0x against a
    1.42–1.57x bound, because Type 1 at chi 0.4 stays one block (its count may not depend on T) and Type 2
    gets floor(T/S) = 1 block per species.
12. A lone `ShiftedCCPotential` with FINUFFT (API use; the `GetCCO` CLI has no `finufft` keyword) now
    costs Type 2 by its largest block, so its gate accepts Type 2 at smaller sizes than before when T > 1;
    forces then differ from earlier builds in the last bits.
13. The `GetCCO` CLI runs Direct sums only, so its contact-cell fix gains about 3%; a `finufft` keyword
    there would give the 1.3–4x evaluation gains measured through the API.

## Provenance

Branches (local, not pushed): `jaeukk/cco-optA-type2-forces` (6e57558, 093a5bb),
`jaeukk/cco-optB-persistent-plans` (5f91e7d, 79b5371, 9f0bb99), `jaeukk/cco-optC-species-concurrent`
(6ce4b8f, 7917475), `jaeukk/cco-optD-indexed-view` (95a608c), and `jaeukk/cco-finufft-integration`
(merges 68b8abf, 4e93bad, a01f36c, 5c11101; fixes 7f5933a to 5d209ca; gate calibration c8147ff). Codex reviews are in
[reviews/](reviews/); task contracts and dispatch records are in `~/.claude/multiagent/tasks/cco-finufft-*`.
