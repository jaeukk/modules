# Results: contact repulsion and fixed-order summation

Completed 72 optimizations: two sum implementations, three cases, three seeds,
two backends, and two independent process repeats. All passed the independent
ground-state checks, with no FINUFFT fallback and successful optimizer pass codes.
See [verification receipt](verification.json) and [per-run timing tables](REPORT.md).

## One component with repulsion

N=500, 299 modes (chi=299/998), packing fraction 0.2, u0=e0=1, two and four
threads. All three initial coordinate arrays and mode files were checked against
the earlier one-component benchmark and match exactly. Contact forces and energy
were validated by independent all-pair calculations and finite differences.

Mean cross-backend Cartesian displacement remained 1.53e-11 to 5.78e-11, including
both baseline and fixed builds. Vanilla repeats were byte-identical. Repulsion alone
therefore did not reproduce the large five-component discrepancy under these conditions.

## Effect of fixed-order summation at four threads

The contact-energy implementation now computes one contribution per particle in
parallel, then sums those contributions serially in particle-index order. The force
formulas, density implementation, and optimizer are unchanged. Baseline executables
retain the original guided OpenMP reduction. No timing runs overlapped.

For N=2000, five components, chi_i=80/798, and packing fraction 0.2:

| Comparison | Original sum: mean displacement range | Fixed sum: mean displacement range |
|---|---:|---:|
| Vanilla repeat versus repeat | 0.01829–0.04357 | 0, byte-identical |
| FINUFFT repeat versus repeat | 0.28139–0.38244 | 0.29247–0.42385 |
| Vanilla versus FINUFFT | 0.27976–0.42359 | 0.29550–0.42899 |

Ranges cover these three seeds, not universal bounds. Earlier repeated runs produced
different displacement magnitudes; that variation is part of the reproducibility issue.
Fixed-build vanilla optimizer pass counts and returned energies also match exactly
between repeats for every seed. This controlled intervention supports the contact
reduction as the source of the observed vanilla nondeterminism in these tests.

In the five-component fixed-point test, vanilla energy range fell from 2.84e-14 to
zero, with zero force variation before and after. FINUFFT energy range was zero after
the change at that sampled configuration, but force variation of 4.44e-16 remained.
The fixed-order contact sum therefore does not make the complete FINUFFT path
reproducible. It also does not explain quantitatively how those residual differences
grow along optimization trajectories; iteration tracing remains unperformed.

All fixed-build five-component energies were between 1.42e-20 and 9.33e-20. All
one-component energies were between 1.69e-22 and 1.27e-21. Every endpoint passed the
same independent energy, force, structure-factor, and overlap checks.

## Median optimization wall time, seconds

Each entry is the median of six runs (three seeds, two repeats). These are full
optimization timings, not identical-path kernel benchmarks; differing evaluation
counts and ordinary timing variation limit before/after speed interpretations.

| Case | Original vanilla | Original FINUFFT | Fixed vanilla | Fixed FINUFFT |
|---|---:|---:|---:|---:|
| One component, two threads | 1.278 | 0.920 | 1.313 | 0.917 |
| One component, four threads | 0.735 | 0.704 | 0.718 | 0.678 |
| Five components, four threads | 6.563 | 10.207 | 7.142 | 10.234 |

The original specialized five-component regression also passed after the change,
including species reassignment/order, cutoff refresh, force finite differences, and
constructed zero-energy configurations at all seven chi masks. `git diff --check`
passed. Code changes are local and uncommitted. The broader sweep remains stopped.

## Remaining interpretation

The data rule out “adding repulsion necessarily causes large endpoint differences.”
They establish that fixed-order contact energy restores vanilla repeatability for
this five-component sample, while FINUFFT repeat variability remains. Lower chi,
component coupling, relative potential strengths, and optimization sensitivity have
not been separated by these tests. No claim of chaos, complete backend equivalence,
or elimination of every source of nondeterminism is justified.
