# CCOptimization: scheduling and L-BFGS history

October 19, 2026. Measurements use implementation commit
`856d6d104dd8a114a63051f61a500899d53e5ac2`; the subsequent default change selects
history 40. [Interactive report](THREADING_RESULTS.html).

## Selected settings

Use **Option A**, **L-BFGS history 40**, and **relative contact strength 100**.
Option A already was master's scheduler and remains the default. FINUFFT remains
opt-in through `finufft`. The shared C++ `RelaxStructure_NLOPT` history default
changes from 100 to 40; `history 100` in GetMultiCCO restores the former setting.
The required contact-strength input must explicitly be 100. The potential still
supports other explicitly requested strengths; no formula or input format changes.

GetCCO and GetMultiCCO run C++ NLopt L-BFGS. Python in the comparison launches
processes and analyzes results; it does not optimize configurations. The old Della
paper adapter separately uses SciPy L-BFGS-B with history 20. Updating its linked
evaluator does not change its optimizer or history. The measurements here do not
establish the performance of that SciPy workflow. New native-C++ sweeps require
an explicit protocol record and separate provenance from the old campaign.

## What the schedules do

**A:** an outer OpenMP pool schedules independent species density work, eligible
force transforms, and contact chunks. Each FINUFFT plan has one internal thread.
The existing particle partition and force cost gate remain in effect. Eligible
Type 2 transforms batch spatial force components. Species count plus one is a
reasonable initial allocation, not a guarantee of full utilization.

**B:** complete each species' density and forces sequentially. Each FINUFFT plan
requests the entire outer thread budget. An additional asynchronous contact
worker time-shares the same allocated CPUs. Four allocated CPUs therefore permit
four internal FINUFFT threads plus one runnable contact worker, not a 3+1 split.
Type 2 transforms batch spatial components without A's concurrency cost gate.
Select this optional schedule with `finufft-sequential`; `finufft` switches back
to A. Compatible plans persist and moving particle positions are refreshed with
`setpts`. Neither schedule requires a different optimizer.

The implementation adds exception-safe contact joining, thread restoration,
direct fallbacks, CLI selection, focused tests, and reproducible comparison tools.
History input now rejects negative, fractional, malformed, and overflowing values
before passing a correction-pair count to NLopt.
The candidate energy in the optional MINOP optimizer is initialized within each
iteration, removing a GCC 11 uninitialized-value diagnostic without changing the
line search or the selected NLopt workflow.
It does not migrate MSPB or SpherePack into core or alter scientific formulas,
the three NLopt scaling phases, acceptance criteria, or storage formats.

## A versus B at contact strength 100

Machine: AMD Ryzen 5 7500F, WSL2, four distinct exposed cores bound to CPU IDs
0/2/4/6. Three equal species, 2D, unit total number density, nominal per-species
chi 0.30, phi 0.70, contact strength 100, FINUFFT 2.5.1 tolerance 1e-12,
history 100. OMP_DYNAMIC=false, OMP_WAIT_POLICY=PASSIVE, OMP_PROC_BIND=false,
OMP_NUM_THREADS=OMP_THREAD_LIMIT=4. Seeds 17–19. Runs are interleaved; timing
excludes setup, warmup, snapshot writing, and post-run direct checks.

Warm evaluation medians use nine batches per schedule and size, 400 fused
`SetConfiguration + EnergyAndAllForce` calls per batch. CPU utilization is total
process CPU time divided by wall time and four allocated CPUs.

| Particles/species | A, ms/call | B, ms/call | B/A wall time | A CPU | B CPU |
|---:|---:|---:|---:|---:|---:|
| 300 | 0.410 | 3.088 | 7.54 | 59.5% | 29.9% |
| 600 | 0.604 | 3.299 | 5.46 | 67.7% | 30.9% |
| 1,000 | 0.805 | 3.558 | 4.42 | 73.8% | 32.4% |
| 4,000 | 2.670 | 6.791 | 2.54 | 83.2% | 41.5% |

Bounded C++ NLopt runs use three seeds per schedule and size, with at most
1,000 objective evaluations per scaling phase. NLopt may overshoot the limit
slightly. These are evaluation limits, not iteration counts.

| Particles/species | A, seconds | B, seconds | B/A wall time | A CPU | B CPU |
|---:|---:|---:|---:|---:|---:|
| 300 | 2.815 | 7.204 | 2.56 | 40.8% | 29.2% |
| 600 | 5.077 | 12.927 | 2.55 | 40.8% | 29.6% |
| 1,000 | 8.093 | 16.213 | 2.00 | 41.2% | 30.1% |
| 4,000 | 28.837 | 35.813 | 1.24 | 38.2% | 31.8% |

All initial and final energy/force checks passed against direct sums. **None of
these 24 bounded optimization runs reached the ground-state threshold.** Their
times are not time to a ground state. Different evaluation counts and trajectories
also limit comparisons of complete optimizations.

A is faster at every measured size. Larger transforms reduce B's disadvantage
here, but these data establish no crossover at which B becomes preferable.
FINUFFT's internal synchronization and small transform sizes are plausible
contributors; no profiler isolates their individual time costs. Serial NLopt
history and line-search work also limit whole-run utilization. Higher CPU
utilization alone is not a performance objective.

An additional 4,000-particles/species bounded replay found byte-identical A
endpoints for the repeated seed 17, but different B endpoints between repetitions.
At identical coordinates, A/B energy and force differences were at most
1.78e-14 and 2.05e-16, respectively. These endpoints were **not ground states**.
Consequently this test cannot establish that B fails to find or reproduce ground
states. Internal floating-point summation order can affect optimizer trajectories;
the specific cause of the observed B variation was not isolated.

## History 40 versus 100 at contact strength 100

Same 4,000-particles/species Option A fixture, phi 0.70, three paired seeds,
three repetitions per seed and history, 1,000 evaluations per phase:

| History | Median wall time | CPU utilization | Median evaluations | Accepted grounds |
|---:|---:|---:|---:|---:|
| 100 | 28.690 s | 38.3% | 2,136 | 0/9 |
| 40 | 14.895 s | 48.8% | 2,035 | 0/9 |

History 40 reduces median wall time by 48.1% (1.93 times faster). Both histories
repeat their own endpoints byte-for-byte for every tested seed. Their endpoints
and final energies differ across histories. Every final scaling phase reported
NLOPT_FAILURE; these timings describe bounded unsuccessful runs. Reduced history
lowers storage and history-update work, but trajectory changes contribute to the
timing difference too. This is not a general ground-state success-rate estimate.

Longer strength-100 pilots at phi 0.50, 0.30, and 0.10 also failed to reach the
acceptance threshold for seed 17 within 10,000 evaluations per phase. No tolerance
was loosened and no failed pilot was classified as an accepted ground state.

## Ground-state identity test: a separate strength-1 experiment

The final short identity test used **contact strength 1**, phi 0.10, Option A,
4,000 particles/species, the same random seed 17, and two runs per history.
It used standard C++ NLopt with 10,000 evaluations per phase and no continuation
or reference-packing starts. The actual finite-system chi was 0.3005751438.

| History | Median wall time | Evaluations | Direct final energy | Accepted grounds |
|---:|---:|---:|---:|---:|
| 100 | 26.997 s | 2,081 | 1.85e-20 | 2/2 |
| 40 | 14.554 s | 2,152 | 1.27e-20 | 2/2 |

All three NLopt phases returned success. Independent direct Fourier sums and
periodic contact checks verified energy <1e-14, constrained species S <1e-10,
maximum force component <1e-7, and relative overlap <1e-6. Across all four
endpoints, maximum energy was 1.85e-20, maximum constrained S was 1.44e-20,
maximum force component was 1.03e-12, and overlap was exactly zero.

Each history repeats its own endpoint byte-for-byte. The two histories produce
distinct valid ground states: labeled periodic RMS displacement after translation
alignment is 1.99 diameters, and translation/permutation/square-symmetry invariant
structure-factor probes also differ. Identical coordinates are not required for
finding degenerate ground states.

The zero-overlap endpoints also pass the same physical checks when independently
reevaluated at strength 100. **They were optimized at strength 1.** This does not
demonstrate strength-100 convergence or a 1.85-times ground-state speedup at
strength 100. Only one random seed was tested. Production remains at strength 100.

## Evidence, reproduction, and verification

Compact raw measurements, summaries, build hashes, environment, and the history
capture harnesses are stored under
[`tests/threading/results-20261019/`](tests/threading/results-20261019/README.md).
Full historical MSPB snapshots and detailed local logs remain in the original
`/tmp/cco-*20261019` artifact directories; they are not committed as research data.
Those local paths are temporary. The committed JSON preserves measured outcomes
and endpoint hashes, but independent coordinate reanalysis requires the snapshots
or a rerun. The report does not silently substitute a new build for the measured one.

The scheduling implementation received two Claude Opus audit cycles before commit
856d6d1. Focused tests cover direct references, finite differences, persistent
plans, motion, batching, thread limits, fallbacks, toggling, and concurrent-contact
exceptions. Release validation also checks the default/explicit history-40 endpoint
identity, history-100 override, CLI export/resume, MSPB archive compatibility, and
offline PairStat. See [threading notes](tests/threading/README.md) for commands.

No Della speedup or full CPU utilization is guaranteed by these local results.
Validate the native execution path and effective settings on a Genoa compute node
before production; keep releases immutable and skip previously explored sweep
points using receipts rather than only accepted-state counts.
