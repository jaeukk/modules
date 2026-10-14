**Independent assessment**

- The remaining evaluator cost is mostly transforms plus contact work; setup is only about 9% of the supplied 2D totals.
- Densities are already reused: `potential/CollectiveCoordinatePotential.h:880` returns cached rho, and `:1019` recomputes only when invalid.
- Contact neighbors are traversed separately for energy and forces (`apps/CCOptimization/MultiStealthyPotential.cpp:383`, `:474`).
- Persistent transforms still execute `setpts` and FFT/spreading work (`potential/CollectiveCoordinatePotential.h:760`, `:845`); removing allocations alone cannot eliminate this floor.
- Full optimization also pays for reconstruction, gradient conversion, and L-BFGS with history 100 (`potential/StructureOptimization.cpp:335`, `:494`, `:869`).
- The evaluator benchmark uses copy-plus-`MoveParticle`, whereas `Objective` uses remove-and-insert reconstruction; it does not measure the exact callback (`tests/bench_finufft.cpp:152` versus `potential/StructureOptimization.cpp:378`).
- The report’s “ranking holds” across release flags is unproven: changing direct-loop performance can change the force-gate crossover (`docs/finufft/perf/REPORT.md:130`).
- Best-of-five timings and loaded-machine medians are different statistics. Neither justifies scaling the entire 15-second optimization by an evaluator speedup.

**Checks performed**

I read the requested report and code, including the cache, force gate, indexed views, optimizer, and build configuration. I ran the existing executables without rebuilding or editing repository files:

| Command, relative to `apps/CCOptimization/` | Result |
|---|---|
| `test_finufft_multi.out` | All FINUFFT multi-species checks passed |
| `bench_finufft.out eval 2 2000 5 0.1 4 1 finufft` | FINUFFT median 1.048 ms; Direct 3.027 ms; relative energy/force differences 4.01×10⁻¹⁵ / 6.39×10⁻¹⁵ |
| `bench_finufft.out eval 2 2000 5 0.4 4 1 finufft` | FINUFFT median 1.308 ms; Direct 9.466 ms; relative energy/force differences 1.48×10⁻¹⁵ / 2.75×10⁻¹⁴ |

Each benchmark invocation internally evaluates both backends for 100 repetitions. These are baseline checks, **not measurements of either proposal**. The disposable copy lacks Git metadata, so I could not independently verify the stated commit or the existing binaries’ build provenance.

**Option 1 — Combine contact energy and force evaluation**

**Mechanism and code.** Add a combined energy-and-force interface to `potential/Potential.h:60`, with a default implementation calling the existing methods. Use it when `Objective` receives a gradient (`potential/StructureOptimization.cpp:455`, `:487`). Override it in `MultiStealthyPotential`, refactoring the separate contact calculations at `apps/CCOptimization/MultiStealthyPotential.cpp:381` and `:470`.

For each particle, traverse its neighbors once and accumulate both its contact-energy slot and its contact force. Keep energy’s `i < SourceParticle` filter and force’s `i != SourceParticle` filter. Reuse the distance calculation while retaining the existing arithmetic expressions and overlap convention.

Keep **directed, per-particle force accumulation**: do not update both particles through atomics. Preserve the ordered energy reductions at `:423`. This removes a neighbor traversal without introducing schedule-dependent sums. The existing FINUFFT transforms, species scheduling, and force gate remain the numerical machinery.

Use the combined path only for explicitly selected FINUFFT evaluation initially. Direct remains the default and fallback. Energy-only requests and cold MD `AllForce` retain their existing entry points.

**Expected gain — estimates, not measured improvements.**

The supplied energy phase includes concurrent species work and contacts, so its full 0.40–0.44 ms is not removable. My planning assumption is that eliminating one contact traversal saves approximately 0.08–0.18 ms at N=2000 in 2D, after allowing for work already hidden behind transforms.

| Case | Estimated saving per gradient evaluation | Evaluator reduction |
|---|---:|---:|
| 2D N=2000, χₛ=0.1 | 0.08–0.18 ms | 8–17% |
| 2D N=2000, χₛ=0.4 | 0.08–0.18 ms | 7–17% |
| 3D N=2000, χₛ=0.1 / 0.4 | 0.10–0.35 ms | 2–9% |
| 2D N=10000, χₛ=0.4 | 0.4–1.2 ms | 2–7% |
| 3D N=100 | 0–0.01 ms | 0–3% |

The 3D and large-N ranges have lower confidence: no contact-only timings were supplied. The large-N estimate allows more than linear scaling of traversal cost, but **does not assume contacts explain the 18.1 ms total**.

For a full optimization, the defensible model is:

`seconds saved = gradient callbacks × milliseconds saved / 1000`

Thus 1000–4000 callbacks imply **0.08–0.72 seconds saved** for either 2D N=2000 case; 3000 callbacks imply 0.24–0.54 seconds, roughly 2–4% of the reported 15-second run. For the other cases, savings per 1000 gradient callbacks are 0.10–0.35 seconds in 3D N=2000, 0.4–1.2 seconds at N=10000, and at most 0.01 seconds at N=100. Their actual full-run savings require their own callback counts; applying the 2D count range to the report’s one-second 3D run would be inconsistent.

**Repeatability, accuracy, and endpoints.** No constraint needs relaxing. Every particle owns its force and energy slot; summation order stays independent of worker assignment. Retaining the original expressions and neighbor order should permit bitwise agreement with the current evaluator, which should also preserve trajectories. This must be tested, not assumed. No FINUFFT tolerance changes are needed.

**Main risks.** A fused implementation can accidentally change periodic-image handling, the coincident-particle rule, or contact counting. Extra force work must not leak into energy-only requests. Restructuring the parallel regions can erase the gain by reducing existing overlap between contacts and species work. Avoid caching across configurations; a combined call needs no new long-lived validity state.

**Validation.**

- Compare combined results against the current separate calls, first bitwise and then against a fresh Direct oracle using the stipulated relative and near-zero bounds.
- Cover 1–3D, skew and small periodic cells with multiple images, unequal radii, coincident particles, disabled repulsion, S=1–10, and FINUFFT failures.
- Repeat across processes at 4 and 8 threads, including passive waiting; exercise energy-only and cold MD force calls.
- Benchmark the combined callback against the separate callback on identical coordinate sequences. Instrument contact traversal separately, then run complete three-pass optimizations and compare callback counts and endpoint hashes.

**Option 2 — Reduce the configurable L-BFGS history**

**Mechanism and code.** Replace the hard-coded history of 100 at `potential/StructureOptimization.cpp:869` with an explicit optimizer setting. Benchmark histories **20, 40, and 100**, selecting a fixed setting offline for production. Thread that setting through `RelaxStructure_NLOPT` and its three passes (`:981`) from the ground-mode configuration in `apps/CCOptimization/main.cpp:1479` and the call at `:1752`.

This attacks optimizer vector work and memory traffic, not potential evaluation. Keep the same three coordinate-rescaling passes, stopping rules, maximum evaluations, and endpoint acceptance criteria. Do not select history from runtime timing.

For particle-only optimization, storing two history-vector families requires approximately `2 × history × dimension × N × 8` bytes, excluding other workspace:

| Case | History 100 | History 20 |
|---|---:|---:|
| 2D N=2000 | 6.4 MB | 1.28 MB |
| 3D N=2000 | 9.6 MB | 1.92 MB |
| 2D N=10000 | 32 MB | 6.4 MB |

History-dependent work can fall substantially once the history fills. This is particularly relevant after accelerating the potential, but it does **not** establish how much of the observed full-run time belongs to L-BFGS.

**Expected gain — explicitly conditional estimates.**

Potential evaluation itself improves by **0 ms in every measured case**. Any improvement is in optimizer overhead per callback and possibly the number of callbacks.

Let `h` be the fraction of baseline optimization time attributable to history-dependent work, and `r` the new-to-old evaluation-count ratio. Assuming comparable iteration-to-evaluation ratios and an 80% reduction of history-dependent work:

`new time / old time ≈ r × (1 − 0.8h)`

The following are planning scenarios to test, not measured overhead fractions:

| Case | Assumed `h` | Full-run change if `r=1.00–1.05` |
|---|---:|---:|
| 2D N=2000, χₛ=0.1 | 5–15% | About 1% slower to 12% faster; −0.04 to +0.47 s saved on 3.9 s |
| 2D N=2000, χₛ=0.4 | 15–30% | 8–24% faster; 1.1–3.6 s saved on 15 s |
| 3D N=2000, either χₛ | 5–15% | About 1% slower to 12% faster; up to 0.12 s on 1 s, or 0.85–1.02 s on 7.1–8.5 s |
| 2D N=10000 | 5–20% | About 1% slower to 16% faster; multiply by measured full-run time |
| 3D N=100 | 0–5% | About 5% slower to 4% faster; little expected value |

For the 15-second case with 3000 callbacks, the assumed history overhead corresponds to approximately **0.6–1.2 ms saved per callback before evaluation-count changes**. That is optimizer overhead, not a reduction of the 1.08 ms evaluator.

The uncertainty is material: with `h=15%`, an evaluation-count increase above about 14% eliminates the predicted gain. I rank this second because it has larger possible full-run value but substantially greater convergence risk.

**Repeatability, accuracy, and endpoints.** A fixed history setting introduces no scheduling-dependent sums and leaves energy/force accuracy unchanged. FINUFFT internal threads remain one; the backend default and fallback are unchanged. Trajectories and endpoint coordinates will change, which the brief permits. Endpoints must still pass the existing ground-state acceptance threshold and independent Direct reevaluation.

**Main risks.** Shorter history may lose useful curvature information, increase line-search work, or cause failures near χₛ=0.4. Early or short runs may never fill 100 history entries and therefore gain little. Reduced optimizer cost is not sufficient evidence: time per successfully accepted ground state is the decision metric.

**Validation.**

- Measure time inside the actual `Objective` callback and total `nlopt_optimize` time separately, with evaluation counts and optimizer status recorded per pass. Do not attribute their entire difference to history operations.
- Compare histories 20, 40, and 100 on identical seeds, initially the four N=2000 cases, then N=10000 and unequal species. Interleave runs on the shared machine.
- Require repeated-process endpoint hashes to match for each fixed setting at 4 and 8 threads.
- Reevaluate every endpoint with a fresh Direct potential; record acceptance rate and time per accepted endpoint across multiple seeds.
- Retain history 100 as the default until the smaller setting wins end to end without degrading endpoint validity or success rate.