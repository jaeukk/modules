# Audit ledger

These are four bounded main-agent test/review cycles, not claims of independent multiagent review. The same seed and independently evaluated full-model acceptance criteria are retained throughout.

## Cycle 1: coordinate scaling — reviewed, no fix

Single L-BFGS runs with coordinate scales 10000, 1, and sqrt(2000) give mean Cartesian cross-backend separations 0.430880, 4.443408, and 0.358722, respectively. The scale-10000 control reproduces the earlier first-pass endpoint behavior and stops above the physical ground-state criteria; it is not the earlier complete three-pass optimization. Both scale-1 and Cartesian-scale pairs pass the independent ground-state criteria, but neither passes the maximum fractional displacement gate of 1e-7. Coordinate scaling materially changes the path but does not fix cross-backend agreement.

## Cycle 2: contact isolation and accumulation — reviewed, no full-model fix

Removing contact repulsion while preserving five components and the initial positions gives mean displacement 1.77986e-12 and maximum fractional displacement 2.34427e-13. This is a CC-only control, not a ground state of the full model: its contact energy under the original objective is about 134.19.

Retaining contact, sorting every particle's contact-force terms by neighbor ID, and accumulating them in long double gives mean displacement 0.362607 and maximum fractional displacement 0.0425298. Both endpoints pass full-model ground-state checks, but their positions disagree. This does not support promoting that experimental change as a repeatability fix. Cycle 4 will independently check the experimental force against the existing implementation, an all-pairs reference, and finite differences.

The initial state has 786 contacting pairs, minimum separation 0.0105981, and pair diameter 0.504627. The contact pair-relative-coordinate radial curvature is +7.85398, while the most negative transverse curvature is -366.112. These are pair-level quantities, not eigenvalues of the whole system. They demonstrate local nonconvex contact geometry but do not measure full-trajectory amplification or prove its exact mechanism.

## Cycle 3: alternative search strategies — reviewed, no fix

Preconditioned restarted truncated Newton reaches its approximately 10000-evaluation limit in both backends without passing ground-state checks; mean separation is 0.745449. It is rejected as a solution in this bounded test. NLopt can slightly exceed its evaluation setting within an inner operation: observed counts are 10014 and 10022.

The second alternative is a diagnostic spectral steepest-descent method in Cartesian coordinates, with Armijo backtracking, bounded Barzilai–Borwein step size, and long-double scalar dot products. The two backends independently supply their own energy and force. No common direct-sum correction, state copying, quantization, or fallback is used. Full contact energy remains enabled. Both runs exhausted 40000 evaluations above ground-state tolerances (energies 1.11e-6 and 1.36e-6). Mean displacement is 0.110246, maximum fractional displacement 0.0255495. This reduces but does not resolve disagreement; no refinement of these already distinct states is promoted as a fix.

## Cycle 4: fixed-parameter search and validation — complete

Since adaptive methods still disagree, the final cycle removes adaptive step selection entirely: fixed Cartesian step coefficient 0.02, velocity retention 0.9, zero initial velocity, and the unchanged full-model force. This is the heavy-ball method. Each backend begins from the original common seed.

After 40000 evaluations per backend, mean Cartesian separation is 3.68697e-12, maximum fractional separation 3.92481e-12. Both energies are about 3.92695e-5, so these are not yet ground states. This is a successful path-agreement control, not a completed solution.

L-BFGS refinement of each backend's own checkpoint produces valid ground states but amplifies separation to mean 6.69974e-5 and maximum fractional 2.36148e-5, failing the 1e-7 position gate. Repeating the refinement exactly reproduces each backend's saved endpoint and energy history. Refinement starts are slightly different, but both derive from exactly the same original seed state; the receipt and summary record this distinction.

The final continuation therefore retains the fixed-parameter method. It restarts each backend from its own 40000-evaluation checkpoint with velocity reset to zero in both runs, permits at most 100000 more evaluations, and stops at energy <=1e-13. The earlier diagnostic stop target was 1e-16 but had not been reached in any fixed/adaptive descent run; the continuation target is still stricter than the unchanged independent ground-state energy criterion of 1e-11. All other physical criteria remain required. Both continuations stop after 84043 evaluations and pass all physical criteria. Mean Cartesian displacement is 3.79989e-12; maximum fractional displacement is 4.02850e-12. This is continuation within cycle 4, not a fifth audit cycle.

An additional short cycle-4 control refines exactly the same late checkpoint with the two backends. It still fails the position gate (mean 7.78171e-5, maximum fractional 2.12367e-5), although both endpoints pass ground-state checks. This excludes slightly different restart coordinates as the sole explanation for late-stage L-BFGS amplification.

Verdict: demonstrated a passing alternative for this original seed, with total optimization times 550.247 seconds vanilla and 501.442 seconds FINUFFT. The original optimizer is not repaired or replaced. Generality across seeds and full end-to-end repeats of the fixed-parameter workflow are not established by this audit. All four cycles are closed; no fifth cycle or broader sweep is launched.

The force validation passed at the original state and saved callback 74: maximum ordered-versus-original force difference 8.88178e-16, maximum difference from the independent all-pairs/direct reference 1.58673e-12, and maximum sampled finite-difference error 1.33834e-9. The finite-difference check covers six particle-coordinate directions per state and backend, not every degree of freedom. Production code remains unchanged.
