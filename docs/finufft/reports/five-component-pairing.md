# Audit of paired ground-state displacements

## Conclusion

The three original vanilla/FINUFFT pairs used byte-identical initial-coordinate files and mode files. The reported large final displacements are real, not a seed mismatch, row permutation, or unit conversion error. However, the four-thread optimizer is not reproducible even with the vanilla backend unchanged. The original comparison did not control for this and cannot attribute endpoint displacements specifically to FINUFFT.

The production sweep remains stopped. No production source changes or fixes were made during this audit.

## Pairing and independent calculation

Every initial file includes all 2000 particle IDs, species, fractional coordinates, and radii. Direct and FINUFFT files match byte for byte within each seed. SHA-256:

| Seed | Initial SHA-256, identical for both backends |
|---|---|
| 20260928 | e0ba16993736d41186a8f61d22686f4063bd3bd7a347b7351fb32e5e954c426c |
| 20260929 | 7535811dcb522b9bb6107d5f0ee20d0549685b572a1cee93427b95eb54f4870d |
| 20260930 | ed98af672f2fdffac5af36ff220cd308424de3f949b15636ff333168f1e931c8 |

All six mode files have SHA-256 e2645ee2b85fcbeac732307913b81011ed40ce894e1cb366070d428f646326f0. The executable also still matches its recorded preparation hash.

Final files each contain 2000 rows in ID order, species=id modulo 5, with matching labels and radii. Independent NumPy evaluation of fixed-label minimum-image displacements, global translation removal, and multiplication by box length sqrt(2000) reproduced the report:

| Seed | Mean | RMS | Maximum |
|---|---:|---:|---:|
| 20260928 | 0.4163172296757778 | 0.5028870519184960 | 2.056904654204547 |
| 20260929 | 0.3951396770761619 | 0.4724827267351626 | 2.1836456757754674 |
| 20260930 | 0.2868330196269472 | 0.34375757382340916 | 1.9468650185154233 |

Global shifts were only about 1e-15 to 8e-15 in fractional units. These values are lengths in the simulation's Cartesian units, not fractional box lengths. There is no particle permutation matching.

## Reproducibility controls

Four isolated repeat optimizations used the unchanged executable, seed 20260928, chi target 0.1, original four-thread environment and CPU affinity 0,2,4,6. Each repeat's initial file matched its original byte for byte. These were diagnostic subprocesses, not a resumed sweep.

| Repeat versus original of same backend | Mean displacement | RMS displacement |
|---|---:|---:|
| Vanilla repeat 1 | 0.36097196357999395 | 0.43728206605251674 |
| Vanilla repeat 2 | 0.35914734784791036 | 0.4362147974532893 |
| FINUFFT repeat 1 | 0.35464686356268404 | 0.4299151654524891 |
| FINUFFT repeat 2 | 0.3768510591699485 | 0.4575751924163632 |

Original vanilla pass evaluation counts were 1051,42,27; the two vanilla repeats used 1026,49,29 and 1037,45,23. Differences therefore appear during optimization, not just in coordinate serialization.

Two vanilla controls with OMP_THREAD_LIMIT=1, OMP_NUM_THREADS=1 and CPU affinity 0 produced byte-identical final files, identical pass counts 1041,61,25, and accepted endpoints. The program still requests four threads internally; the runtime limit serializes execution.

A FINUFFT control under that same forced runtime limit crashed with SIGSEGV. No FINUFFT single-thread endpoint comparison is available. This thread-limit/request mismatch requires a separate investigation or a properly configured single-thread harness; it is not evidence that the original four-thread runs crashed. The incomplete serial_finufft_0 directory is retained.

## Evidence for nondeterministic energy evaluation

A diagnostic-only driver, fixed_point.cpp, includes the original experiment and rebuilds potential configuration/caches at the same coordinates 100 times, evaluating energy and forces using four threads:

| Backend | Minimum energy | Maximum energy | Energy range | Maximum force change |
|---|---:|---:|---:|---:|
| Vanilla | 137.22535335512228 | 137.22535335512231 | 2.8421709430404007e-14 | 0 |
| FINUFFT | 137.22535335512271 | 137.22535335512273 | 2.8421709430404007e-14 | 4.4408920985006262e-16 |

MultiStealthyPotential.cpp:243 uses an OpenMP guided-schedule reduction for contact energy. Floating-point summation order is not fixed there. The collective-coordinate energy instead uses static thread partitions followed by an ordered sum of thread partials. The NLopt objective consumes both energy and forces; identical gradients alone do not guarantee identical optimizer decisions when the energies differ.

The contact-energy reduction is a concrete source of nondeterministic arithmetic and a strong candidate for the observed trajectory divergence. The audit did not record iteration-by-iteration trajectories or isolate that reduction with a deterministic replacement, so it does not prove this is the sole cause or locate the first divergent line-search decision. FINUFFT also has approximation differences independent of threading. Same initial coordinates plus a deterministic optimizer do not imply identical trajectories across nonidentical numerical objective implementations.

## Interpretation and next corrective test

The existing independent endpoint checks support low-energy, low-force, nonoverlapping configurations; they do not establish deterministic path equivalence. The original endpoint displacement table is numerically correct but is confounded as a backend comparison. Timing measurements are actual elapsed times, but differing trajectories/evaluation counts mean they are not pure kernel speed comparisons.

Before drawing backend-equivalence conclusions, use a fixed-order contact-energy sum, verify same-backend four-thread repeatability, then compare paired iteration traces and FINUFFT tolerances. Those production changes were not made because this request was an audit, not authorization to implement a fix.

Diagnostic coordinates and fixed-point driver/binary are preserved in pair-audit-artifacts/. Original results were not overwritten. The original report's running status is stale because its process was stopped externally; systemd verifies inactive with MainPID=0.
