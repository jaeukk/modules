# Independent validation discovery — September 22, 2026

Status: bounded discovery complete; Gate A and Gate B remain **NOT RUN**. No backend implementation, optimization, full experiment, package installation, or source modification was performed. Owned outputs are confined to `validation/`. Read the entire launch BRIEF and `control/COORDINATOR.md`; the latter is under the task directory, not the workspace root.

## Evidence and provenance

Inspected original source read-only at `/home/jaeukk/30_Codes/modules`, HEAD `ccdd9285b5f5cac37d37766bfa37f98cfd84ee3c`. The source worktree has unrelated modifications and untracked files; HEAD alone is not sufficient provenance. `direct_probe_results.json` records SHA-256 of the five inspected numerical source files. The repository worker owns the full committed-history clone at `../source/modules`. At the final checkpoint, SHA-256 comparison confirmed all five inspected numerical files exactly match that clone; the formulas and protocol therefore apply to the clone. The probe now reads provenance from the clone. Unrelated dirty original files were not applied by this worker.

Primary source locations (line numbers in the inspected original):

- `potential/CollectiveCoordinatePotential.h:82–201`: positive Fourier sign, direct sum, energy, forces, cache.
- `apps/CCOptimization/MultiStealthyPotential.cpp:17–45, 177–332`: species mapping, subconfiguration construction, weights, contact energy and force, parallel wrapper.
- `apps/CCOptimization/main.cpp:1055–1241`: random packing, radius scaling, species assignment, chi input and actual reporting; `1373–1433`: default acceptance energy and limits; `1651–1698`: minimizer dispatch and final success classification; `1856–1864`: continuum chi conversion.
- `core/PeriodicCellList.cpp:509–635`: reciprocal selection and half-space convention.
- `potential/StructureOptimization.cpp:850–898`: NLopt LBFGS storage 100, max evaluations, stop value, status handling.
- Project `30_Codes/pixels/avian_photoreceptor.py:39–59` and `20_Data/color_mixing/2026-09-21-avian-photoreceptor/cco_validation.json`: existing five-component phi=0.20 evidence.

Executable `direct_probe.py` uses only Python's standard library. It independently checks the formulas below; it does **not** invoke compiled C++ or FINUFFT. Reproduce from the workspace root:

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 timeout 30s python3 _multiagent/tasks/2026-09-22-graft-finufft/validation/direct_probe.py
```

Recorded results: best central-difference absolute error 1.61e-11 (h=1e-5), wrapped rho error 8.90e-15, translated energy error 1.11e-16, translated gradient error 1.94e-16, full ±k/half-space energy ratio 2.0000000000000004, and a 4×4 cancellation lattice energy 9.60e-31. All probe assertions passed. This is derivation evidence only, not evidence that Gate A passed.

## Independent formulas and interpretation

Let A have physical cell vectors as columns, V=|det A|, fractional position t_j, physical position r_j=A t_j, and integer reciprocal label m. Then k_m=2π A⁻ᵀm, k_m·r_j=2π m·t_j. For FINUFFT coordinates use x_j=2π wrap(t_j) in a documented length-2π interval; subtracting 2π times integers preserves phases. Simply shifting [0,1) coordinates by 1/2 introduces a mode-dependent phase and must not be confused with wrapping.

For species s with N_s particles, define

- rho_s(m)=Σ_{j∈s} exp(+i k_m·r_j)=R_s+i I_s.
- S_s(m)=|rho_s(m)|²/N_s.
- E_s=V⁻¹ Σ_{m∈H_s} w_sm |rho_s(m)|² + Shift_s.
- E_CC=u0 Σ_s (N_s/N) E_s.
- ∂E_CC/∂r_jα = (2 u0 N_s/(N V)) Σ_m w_sm k_mα [I_s cos(k_m·r_j) − R_s sin(k_m·r_j)].
- Force is minus this gradient. For fractional optimization variables, ∇_t E=Aᵀ∇_r E; the optimizer's additional coordinate rescaling must also be applied.

The executable source, rather than the introductory header comment, establishes **1/V**, not N/(2V), as the base energy factor. H_s contains one representative of each ±k pair. Retaining that half-space gives the formula above; populating both halves at full weight doubles energy and gradient. An adjoint/type-2 evaluation of coefficients w_sm k_mα rho_s(m), with the appropriate opposite Fourier phase, supplies the reciprocal-to-particle force sum; it is not inverse reconstruction of positions.

For the 2D study with diameter below half the box length, the contact term is

E_CP=e0 Σ_{i<j} max(0, 1−d_ij/(a_i+a_j))²,

where d_ij is periodic separation. For 0<d_ij<sigma_ij, with delta=r_i−r_j in the selected periodic image, its gradient contribution on i is −2 e0 (1−d_ij/sigma_ij) delta/(sigma_ij d_ij); j gets the opposite. The actual neighbor iterator sums eligible periodic images, so the minimum-image formula must not be extrapolated to oversized radii or tiny cells. At exact coincidence the radial derivative is undefined; the existing force selects an arbitrary x direction. Test finite behavior there, not equality to a nonexistent classical derivative.

Five components mean five disjoint species of 100 particles each, independently constrained Fourier densities, plus contact interactions among **all** particles (both same-species and cross-species). Each species Fourier energy has weight 1/5. There is no explicit total-density Fourier objective or explicit cross-spectrum penalty. On a common constrained mode, rho_total=Σ_s rho_s; making every rho_s zero therefore makes the total and all cross terms zero automatically. A total-density-only zero does not enforce species stealthiness. If cutoffs differ, this implication holds on their common intersection only. Species count and radius-group count are separate inputs; five species need not mean five distinct radii.

Use S0=0 and flat positive weights for these gates. The multi-species driver accepts S0≠0 but constructs `ShiftedCCPotential`; its `AddConstraint` uses only vals[0], so the intended S0 target is ignored in this path. Do not represent this branch as a validated nonzero-S target.

## Mode counts and chi

GetKs selects the first nonzero reciprocal integer component positive, excludes zero, and sorts by squared physical norm. In 2D the half-space is m_x>0 or (m_x=0 and m_y>0). The driver's lower cutoff is strict. For boundary-sensitive modes explicitly compare the neighbor-iterator and enumeration paths; avoid numerical ambiguity by putting study cutoffs strictly between shells and save the actual k list in full precision.

Report chi_s=M_s/[d(N_s−1)], with M_s the retained nonzero mode-pair count, not the number of real equations. Each complex condition normally supplies two real equations; degeneracies can reduce their rank, so chi is a conventional count, not a proof of independence. For the CC-only product of species, translations provide d zero modes per species; cross-species contact generally leaves only the common d translational invariances. Record both Σ M_s/[d(N−C)] and Σ M_s/[d(N−1)] if discussing aggregate constraint fractions, and name the denominator; neither replaces chi_s.

The driver input uses total number density rho=N/V through K=2π[2d rho chi_input/v_d]^(1/d). It is a continuum cutoff parameter, not realized finite-size chi_s. Approximate conversion to desired per-species chi is chi_input≈chi_s(N_s−1)/N; passing 0.3 directly for each of five species gives approximately chi_s=1.5, not 0.3.

For the square cell V=N=500, use complete shells nearest the requested chi and a cutoff at the midpoint in squared integer norm between shells. Independent integer enumeration predicts:

| Case | Target chi_s | M_s | Realized chi_s | Midpoint q²=(K L/2π)² | Driver chi_input |
|---|---:|---:|---:|---:|---:|
| One species, 500 | 0.30 | 300 | 0.3006012024 | 193.5 | 0.3039490892 |
| Five species, 100 each | 0.10 | 18 | 0.0909090909 | 11.5 | 0.0180641578 |
| Five species, 100 each | 0.20 | 40 | 0.2020202020 | 25.5 | 0.0400553063 |
| Five species, 100 each | 0.30 | 60 | 0.3030303030 | 38.5 | 0.0604756586 |
| Five species, 100 each | 0.40 | 80 | 0.4040404040 | 51.0 | 0.0801106127 |
| Five species, 100 each | 0.50 | 98 | 0.4949494949 | 64.5 | 0.1013163631 |

The 0.10 target ties M=18 and M=22 in exact distance; select the lower shell by convention. Record these as predictions until the actual application exports and confirms its modes. Splitting shells to hit a requested number exactly loses square symmetry. The literal single-component input 0.30 instead yields M=296 and chi=0.2965931864; it is also an acceptable approximately-0.3 benchmark if explicitly chosen and matched across backends.

## Packing fraction decision

For the random driver V=N, giving number density 1. Radii are rescaled so phi=Σ_s N_s v_d(a_s)/V. Equal 2D radii are a=√(phi/π); for arbitrary V they are √(phi V/(Nπ)). Positive e0 makes radii soft exclusion distances, not a hard non-overlap guarantee. At zero contact energy every relevant pair is non-overlapping; at positive energy overlap is allowed. With e0 disabled (|e0|≤1e-10), radii only decorate the point configuration and phi is nominal disk area, not the union covered area if disks overlap.

**Recommend phi=0.20 as the fixed five-component study value**, following the existing five-component precedent rather than falling back to exploratory 0.30. The September 21 validation artifact reports five N=3600 states with 720 particles per species, radius 0.2523132522, 288 pairs per species, actual chi_s≈0.200278, and minimum separations close to twice that radius. This is inspected prior evidence, not an independent rerun or proof for N=500 or higher chi. Other project scripts use phi=0.50 for different component mixtures; there is no universal project phi.

Exploratory phi=0.30 is mathematically coherent (a≈0.309019, nominal area sum 0.30), but is not necessary as a missing-definition fallback and could make joint constraints harder. Use e0=u0=1 as an explicitly documented protocol choice unless baseline configuration establishes another positive ratio. The principal transform benchmark can use e0=0 to isolate Fourier performance, provided both backends use the same value; then add e0=1 wrapper coverage before Gate B and the five-component branch.

## Gate A: required acceptance protocol, not yet executed

Use float64, a fixed cell, fixed species and exact same exported modes/weights. Test eps=1e-6, 1e-9, 1e-12, and 1e-14; only the tight validated setting may drive ground-state experiments. FINUFFT's requested eps is not a formal per-mode error bound, so use the following explicit empirical contract and investigate any breach, rather than assuming the request guarantees success.

For each species set eta=max(50 eps, 200 u N_s), where u=2.22e-16, and delta_s=N_s eta. Require max_m |rho_fast−rho_direct|/N_s≤eta and RMS(|delta rho|)/N_s≤eta; also record relative L2 error when the direct norm is nonzero. Empty mode sets return zero Fourier energy and forces without plans. Evaluate near-cancellation cases by absolute error, not relative error divided by nearly zero rho.

Let a_s=|u0|N_s/(NV). For a direct state define

B_E=Σ_s a_s Σ_m |w_sm| (2 |rho_direct,sm| delta_s + delta_s²).

Require |E_fast−E_direct|≤B_E plus 200u times the sum of absolute direct energy contributions (including any shift and contact terms), with a tiny absolute numerical floor 1e-28 in the chosen energy units. Separately check the energy assembled from fast rho against its reported value, preventing a compensating normalization error.

For each species define Q_s=Σ_m |w_sm| |k_m| (|rho_direct,sm|+delta_s), D_s=Σ_m |w_sm| |k_m| delta_s. A conservative gradient component budget is B_g,s=2a_s(D_s+eta Q_s), plus 200u times the direct absolute summand scale. Require every force component to agree within this bound and report normalized L∞ and RMS errors. These budgets reflect forward density and adjoint error; do not use a fixed relative gradient tolerance near a stationary point.

At the production eps, additionally require B_E<0.1 E_accept and B_g,s<0.1 g_accept on cancellation/near-ground-state fixtures. Proposed fixed benchmark criteria are E_accept=1e-14 for nonnegative shift-free total energy and g_accept=1e-8 in physical coordinates with density 1 and u0=e0=1. If errors cannot be separated from these thresholds, mark the case numerically inconclusive; tighten eps or investigate, and do not silently relax acceptance. Track energy contributions separately to prevent shift cancellation from creating false ground states.

Required matrix:

1. Tiny N=1,2,7,16; N=500 nonoptimized fixture; one and five species, unequal populations and nonuniform positive weights; rectangular and skew cells. Compare exported rho, energies, physical forces, fractional gradients, and component-to-global scatter. Test signed particle weights only if the selected backend exposes that existing weighted class.
2. Random states and analytic cancellation lattices, near-cancellation perturbations, clustered particles, zero/one modes, empty species, zero weights, manually inserted zero mode, duplicate and conjugate modes, isolated high integer modes, and odd/even grid dimensions with boundary index modes. Unsupported cases must give an explicit error or documented direct fallback, never silent wrong answers. The zero mode has rho=N_s and zero gradient; it must not accidentally enter the stealthy mask.
3. Independent per-particle integer box shifts, negative/out-of-box fractional coordinates, points just on either side of a periodic face, common translations, allowed cell symmetries, and within-species permutations. For translations compare rho after its expected phase rotation; energies and corresponding forces are invariant. Noninteger arbitrary k does not obey periodic wrapping: reject or use a compatible fallback.
4. Central differences for every coordinate of N≤16 and seeded directional checks at N=500: h×max(1, typical spacing) with h=1e-3,1e-4,1e-5,1e-6. Require an O(h²) trend until roundoff/transform noise and a minimum error ≤1e-6 times a nonvanishing gradient scale, plus 10 times the predicted evaluation noise divided by h. Use tight eps and ensure the noise term is materially below the criterion. Validate direct gradients first; test random smooth states, avoiding exact overlap/contact transitions. Near stationary points use absolute analytic fixtures instead.
5. Repulsion alone, CC alone, and their sum; pairs below/at/above contact and across periodic boundaries; exact overlap finite-value handling separately. Verify total-force cancellation. Changing radius must update the neighbor cutoff before contact checks.
6. Repeated Energy→Force and Force→Energy, cold AllForce before Energy, SetConfiguration after a move, changed constraints, changed species mapping, repeated plan reuse, and one versus two threads. Test cache invalidation and concurrency, not just the transform kernel.

Read-only code risks requiring focused tests, not asserted runtime failures: `MultiStealthyPotential::AllForce` can enter inner Force/GetRho concurrently on cold caches; species subconfigurations are built in global-index order while stored within-species indices follow assignment order, so arbitrary reassignment needs verification; base AddConstraint does not visibly invalidate cached rho; SetSecureRadii does not update rc. The normal sequential initial species assignment avoids the ordering issue, but unsupported mutations need a clear contract. Any preexisting baseline failure must be distinguished from FINUFFT regressions.

## Matched benchmark and Gate B

Primary: N=500, d=2, square L=√500, one species, flat weights, zero shift/S0, M=300 above (or explicitly choose the literal-input alternative), float64. Primary seed 20260922; paired additional seeds 20260923–20260926. Generate and save initial coordinates once per seed with generator/version and SHA-256; read the same artifact for both backends. The driver's setup consumes an initial random configuration while defining the potential, so matching a seed alone is insufficient. Save exact modes, initial arrays, radii, labels, cell, backend revision, build flags, libraries, eps, thread/plan options, and command.

Use the same existing fixed-cell LBFGS, vector storage 100, maximum 10000 evaluations, and fmin=0 passed by the current caller, unless both arms deliberately receive an identical documented change. The input `tolerance` defaults to 1e-14 and classifies the final result; it is **not** the NLopt stopping value. Max-evaluation termination or an NLopt error cannot be labeled optimizer convergence merely because a Success file exists. Capture raw status and evaluation count explicitly; distinguish achieved ground-state quality from optimizer termination reason. Match thread count (one and two as separate comparisons), precision, initial state, and all stopping/limit settings.

Always direct-recompute final energies, rho/S, forces, maximum and RMS |rho|/N_s, mean and maximum S_s, separate contact energy, and worst relative overlap. Ground-state quality requires direct total energy <1e-14, physical max gradient <1e-8, and separately bounded contact/CC contributions. For positive flat weights each constrained |rho_s|/N_s must also satisfy √(E_accept/(a_s N_s²)); report the actual maximum and aggregate, not just this loose bound. Report the sum of all constrained S_s and their maximum. With unequal or vanishing weights apply unweighted residual criteria too.

Compare periodic positions using minimum-image displacement after an optimized common translation; allow square D4 rotations/reflections only when the exact mask/weights preserve them, and permutations only within identical species. Compare label-preserving and assignment-matched RMS/max distances and retain unaligned coordinates. Independent species translations are symmetries only when contact is absent. Different positions can represent degenerate ground states: compare constrained structure factors, total and partial S(k) outside the mask, and total/partial radial pair statistics with identical bins and finite-size conventions. Five seeds describe variability; do not claim equality of ensembles from that sample size.

Timing: no coincident experimental job during measured benchmark runs. Record monotonic wall time and process CPU time for setup/plan, set-points, forward transform, adjoint/force, uncached energy+gradient, cached repeat, contact, and entire optimization. Count all evaluations and transforms. Warm up each backend with three untimed evaluations, then collect seven timing batches, each at least 0.1 s where practical and capped at 5 s. Alternate arm order across seeds; report medians and ranges, cold and warm costs, and direct/FINUFFT speedup without assuming improvement. Measure peak RSS via `/usr/bin/time -v` in separate matched processes, plus retained plan memory if exposed. Include serialization/setup in a separate end-to-end number. Optional post-gate size sweep N=100,500,2000 at comparable realized chi is capped at 15 minutes total and omitted if resources do not justify it.

Gate B requires Gate A evidence plus all five paired seeds classified, no unexplained backend-specific failures, and comparable direct-checked ground-state quality. A negative NLopt status with a verified candidate should be labeled explicitly and investigated before Gate B, rather than hidden. If both arms miss the acceptance target, report unsuccessful effort and defer the dependent study; do not infer infeasibility or silently count matching failures as comparable converged quality.

## Experiment and resource budget

Observed 2026-09-22 approximately 18:03 Asia/Seoul: WSL/Microsoft VM on AMD Ryzen 5 7500F, 8 visible logical CPUs (4 visible cores, 2 threads/core), 28 GiB RAM with about 26 GiB available, 8 GiB unused swap, 135 GiB free disk, load averages 0.89/0.34/0.11. This is a snapshot, not a reservation. The independent probe finished within its 30-second timeout; no optimization runtime estimate is available yet.

Every computational process: OMP_NUM_THREADS=2, OMP_THREAD_LIMIT=2, OMP_MAX_ACTIVE_LEVELS=1, OPENBLAS_NUM_THREADS=2, MKL_NUM_THREADS=2, NUMEXPR_NUM_THREADS=2, and explicit application/FINUFFT threads≤2. Avoid nested thread pools. At most two concurrent experiment jobs; time benchmark arms serially. Reinspect available memory and load before launch; postpone large cases if available memory <8 GiB or load indicates contention. Suggested per-job RSS budget 4 GiB; record actual maxima rather than relying on input size alone.

After implementation, allow Gate A jobs 60 seconds each and 15 minutes total; initial paired optimization pilot 10 minutes per arm, 10000 evaluations. Use measured pilot cost to forecast remaining paired seeds and 25 five-component starts before launching. Initial study cap 30 minutes per seed and 12.5 wall-job-hours total; if cost is higher, coordinator decides a smaller documented completion set rather than launching unbounded runs. These are budget proposals, not completed experiments or promised runtimes. Use `timeout --signal=TERM --kill-after=30s` plus application checkpoints; a hard timeout is not convergence and may prevent final-state saving.

After Gate B: five independent random-start seeds 20260922–20260926 at each common chi target, N_s=100, phi=0.20 as above. Reusing seed identifiers across chi permits paired initial configurations but preserves independence within each point. Use an independent seed stream for additional starts. Refinement/continuation is a separately labeled branch with parent-state hashes; do not pool it with random starts. Direct-check every N=500 final candidate where inexpensive, at minimum all claimed successes and representative failures. Outcomes: candidate found, not found within allotted effort, or numerically inconclusive. Neither a failed minimization nor no success among five starts proves infeasibility.

Collection cutoff: 2026-09-27 12:00 Asia/Seoul (03:00 UTC). Hard deadline: 18:00 Asia/Seoul (09:00 UTC). Before every job set its timeout to the smaller of its individual cap and time remaining to collection minus a 10-minute save margin; do not launch if this cannot accommodate the forecast. Reserve the final six hours for collection/reporting. The coordinator owns continuation and global concurrency; this discovery worker settles after its report.

## Remaining work and blockers

No scientific-definition blocker remains for the equal-radius five-species S0=0 study using established phi=0.20. Implementation, full-clone baseline, version-matched FINUFFT/Graft evidence, compiled-source comparisons, actual mode export, production-eps selection, and both gates are outstanding dependencies. Source-level cache/mapping/contact-cutoff risks above require baseline tests or scoped fixes by the assigned implementation owner. No Gate A, Gate B, speedup, N=500 ground state, or five-component feasibility result is claimed here.
