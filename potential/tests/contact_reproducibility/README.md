# Contact-energy reproducibility

This finite regression compares optional FINUFFT densities with direct sums in the
presence of contact repulsion. It does not restart the scientific sweep.

Cases: N=500, one component, 299 modes, two or four threads; N=2000, five
components, 80 modes per component, four threads. Total packing fraction is 0.2,
number density is 1, and both potential coefficients are 1. Each case uses three
seeds and two repeats per backend. Requested and allowed thread counts match.

`experiment.cpp` reuses the five-component harness with compile-time test sizes.
The ordinary five-component build retains its original defaults. The original
`validate` branch is a specialized N=2000/five-species fixture, not a general-size
validation; use this driver's `probe` for the other sizes.

`probe` compares energy and forces to an independent direct/all-pair reference,
checks finite differences, and records ranges over 100 repeated fixed-coordinate
evaluations. The runner saves coordinates, independent endpoint checks, optimizer
codes/counts, backend-use counts, input hashes, and optimization timings. It asserts
identical starting coordinates and modes before comparing endpoints. The one-species
generator and 299-mode ordering match the earlier one-component benchmark.

Run `python run.py baseline --output NEW_EVIDENCE_DIR --reference-build ORIGINAL_BUILD_DIR`
before changing the contact-energy reduction, then run the same command with `fixed`
afterward. The reference build must be the five-component build with `commands.json`
and `source_hashes.json`; shared source hashes are checked before reusing its core
objects. Both stages recompile the current contact potential and test drivers.
Never label a build made after the fix as the old baseline. Each phase refuses to
overwrite an existing directory. Run phases serially for meaningful timing.

After both phases, run `python verify.py EVIDENCE_DIR`. This checks all 72 endpoint
receipts, optimizer codes, actual backend use, input identity across phases, fixed-point
vanilla repeatability, and byte-identical fixed-build vanilla endpoints/pass receipts.
It saves `verification.json` and exits nonzero on any failed check. FINUFFT endpoint
identity is measured, not asserted: the transform may introduce its own roundoff variation.

The fixed implementation stores each particle's pair-energy contribution separately,
then sums in particle order. Merely selecting a static schedule with an OpenMP
reduction would not fix the combination order. This change stabilizes the contact sum;
it does not promise deterministic FINUFFT transforms or identical arithmetic across
different thread counts in other potential terms. Force formulas are unchanged.

## FINUFFT-only thread control

`run_finufft_threads.py --fixed-results EVIDENCE_DIR/fixed --output NEW_DIR`
relinks the preserved five-component fixed-build objects with the diagnostic
`finufft_threads.cpp` wrapper. The wrapper copies the supplied FINUFFT options and
changes only `opts.nthreads` to 1. It covers both `finufft2d1` (the preserved objects) and the
persistent-plan path of current builds (`finufft_makeplan`, then one `finufft_execute` per
transform). All other OpenMP work and the runtime thread limit
remain at four threads. Production defaults and the shared library are not changed.

The driver performs 100 fixed-point evaluations per backend, finite differences,
and twelve optimizations (three seeds, two repeats, two backends). Every actual
FINUFFT call is counted, with checks for the expected incoming four-thread request,
the explicit one-thread override, and successful return statuses. Inputs must match
the preceding fixed-build experiments byte for byte. The output report distinguishes
same-backend repeatability from cross-backend endpoint agreement and includes timings.
