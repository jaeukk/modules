# Identification of the requested Zhang comparison protocol

The three requested one-seed diagnostics are complete; see REPORT.md. Cross-backend endpoint agreement remains unresolved.

The local read-only Zotero API identifies citation key `zhang_perfect_2016` as item `BLJC3BSF`: G. Zhang, F. H. Stillinger, and S. Torquato, *The Perfect Glass Paradigm: Disordered Hyperuniform Glasses Down to Absolute Zero*, Scientific Reports 6, 36963 (2016), DOI 10.1038/srep36963. The linked PDF attachment is `9BXFHMGV`.

## Primary-source findings

The [main article](https://www.nature.com/articles/srep36963), Methods, uses low-storage BFGS followed by MINOP, reporting a final gradient norm below 1e-13. Its Discussion checks whether small-system candidate ground states are related by translation, rotation, or inversion. Methods also compares structural statistics for different cell shapes. The [supplement, section 5](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fsrep36963/MediaObjects/41598_2016_BFsrep36963_MOESM1_ESM.pdf) checks agreement between a translational-order metric calculated from pair correlations and from structure factors.

These are distinct procedures. The user subsequently clarified that the intended criterion is geometric identity of configurations, not the minimization procedure.

## Scope and status

- A geometric-identity rejection test has now been performed, as detailed below. No MINOP refinement or structural-statistics comparison was performed.
- The paper uses a different physical potential; applying its numerical checks must not silently replace the current five-component stealthy/contact objective.
- Structural-statistics agreement alone would not establish identical particle configurations or resolve cross-backend trajectory repeatability.
- No Zotero library entries or notes were changed.

## Geometric identity test, September 28, 2026

We apply the paper's geometric-equivalence idea to the saved vanilla and FINUFFT endpoints from seed 20260928. Input arrays confirm exactly identical initial states. Rather than rely on a local alignment search, we use a necessary invariant: the complete sorted list of minimum-image pair distances. Translation, rigid rotation/inversion of the periodic configuration and cell, and particle permutations preserve this list. For a fixed square cell, allowed orthogonal lattice symmetries comprise its eight signed coordinate permutations. Rotating particles arbitrarily while holding the lattice fixed is not generally a symmetry of this periodic problem.

This distance-spectrum test is our implementation of a necessary condition, not an algorithm or numerical threshold attributed to the paper. If an aligned particle matching has maximum positional error e, the triangle inequality bounds each pair-distance error by 2e. Sorting cannot increase this bottleneck bound. Therefore half the largest sorted-distance gap is a lower bound on the best possible maximum positional error. Spectrum agreement alone would not prove geometric identity; a failed necessary condition is sufficient to reject it.

The existing audit tolerance is 1e-7 in fractional Euclidean position units. All distances below use the box side sqrt(2000) for Cartesian conversion.

| Selection | Pair count | Largest sorted-distance gap, fractional | Lower bound on best maximum position error, Cartesian |
| --- | ---: | ---: | ---: |
| All particles, ignoring component labels | 1,999,000 | 0.000421545 | 0.00942603 |
| Component 0 | 79,800 | 0.001699595 | 0.03800410 |
| Component 1 | 79,800 | 0.002022489 | 0.04522424 |
| Component 2 | 79,800 | 0.002297465 | 0.05137287 |
| Component 3 | 79,800 | 0.002476498 | 0.05537617 |
| Component 4 | 79,800 | 0.001789331 | 0.04001066 |

All six tests reject identity. Even discarding component labels and allowing arbitrary relabeling, the lower bound on maximum fractional position error is 0.000210772, over 2100 times the audit tolerance. Each component also differs internally, even when considered independently under rigid transformations and within-component relabeling. These are lower bounds on maximum error, not measured optimal RMS or mean displacements; they do not replace the earlier fixed-ID mean of 0.430886.

Validation: eight synthetic square-lattice orthogonal transformations, each combined with translation, wrapping, and random particle permutation, preserve spectra to better than 1e-14. A deliberately displaced-particle control is rejected. The test reads saved endpoints only and changes neither production code nor optimization results.

Conclusion: these two endpoints are geometrically distinct periodic configurations, not merely different coordinate representations or particle labelings. This resolves the identity question for this pair, not the underlying cross-backend trajectory sensitivity.

Reproduction: `python _multiagent/tasks/finufft-energy-diagnosis/identity.py`. The script refuses to overwrite its existing output. [Machine-readable results and input hashes](results/geometric_identity.json).
