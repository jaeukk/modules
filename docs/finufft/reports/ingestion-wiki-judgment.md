# FINUFFT wiki-ingestion judgment

## Executive decision

Adopt **12 atomic wiki tasks**: five new pages and seven additive updates. This is ready for one-high-tier-worker-per-page dispatch; the exact claim/source contract is in `page-claim-source-map.json`, and the retained duplicate evidence is in `duplicate-find-receipts.json`.

The five new owners are:

1. `20_Topics/01_Math/12_LinearAlgebra/Nonuniform_Fourier_Transform.md`
2. `20_Topics/01_Math/12_LinearAlgebra/Exponential_Semicircle_Kernel.md`
3. `20_Topics/01_Math/12_LinearAlgebra/NUFFT_Aliasing_Error.md`
4. `20_Topics/01_Math/12_LinearAlgebra/Parallel_Spreading.md`
5. `20_Topics/00_Entities/FINUFFT.md`

The seven existing owners to extend, without rewriting or removing anything, are:

1. `20_Topics/01_Math/12_LinearAlgebra/Fast_Fourier_Transform.md`
2. `20_Topics/01_Math/12_LinearAlgebra/Least_Squares.md`
3. `20_Topics/01_Math/12_LinearAlgebra/Inner_Product_Space.md`
4. `20_Topics/01_Math/04_Calculus/Numerical_Quadrature.md`
5. `20_Topics/12_Waves_and_Acoustics/Harmonics_and_Fourier_Series.md`
6. `20_Topics/10_Classical_Electrodynamics/05_Electrostatics/Poisson_Laplace_Equations.md`
7. `20_Topics/02_Statistical_Physics/10_Models/Gaussian_Random_Field_Media.md`

Both the paper and manual source gates are closed. The repaired manual **revision 2** passed full independent re-audit (`numerical/ingestion/manual-reaudit/REPORT.md`, SHA-256 `76f595627254b03447d53ecfc199889a9915f99cbd067a646521330369227ff4`); `RAW_BUILD_COMPLETION.md` and `acceptance/manual-acceptance.json` confirm all 104 files were published and read back canonically with zero broken manual links, zero soft links, zero alias-target collisions, and zero new trunk breakage. I rechecked every cited manual block ID in the canonical notes. There is therefore no remaining source-audit dependency blocking page dispatch or publication; only the ordinary per-page implementation and conductor-owned final wiki checks remain.

## Source and read boundary

The vault was resolved at runtime from `~/.claude/zotero-obsidian.local` as `/mnt/c/Users/김재욱/My Drive/_WORKSPACE/20_Notes`. I read the controlling vault instructions, policy, wiki schema, math notation, and actual wiki-ingest skill, plus the task brief, numerical handoff/report, manual route/destination, build provenance, coverage tree, manifest, acceptance receipts, and current independent audit.

The dispatched worker configuration was checked through Orca: requested and effective model `gpt-5.6-sol`, reasoning effort `ultra`.

### Accepted papers

- `40_Resources/10_Papers/Simulations/barnett_parallel_2019.md`: read completely. This is the accepted published SIAM paper, DOI `10.1137/18M120885X`, pp. C479–C504, Zotero parent `XBQ7PTEV`, attachment `VJBBMK2G`. The original verifier's exit 2 is a layout abstention; the coordinator-approved independently rendered source-label adapter passes. The complete source note has 42 numbered equations, 219 unnumbered math occurrences, and two full tables.
- `40_Resources/10_Papers/Simulations/barnett_aliasing_2021.md`: read completely. This is the exact published ACHA paper, DOI `10.1016/j.acha.2020.10.002`, 16 pages, Zotero parent `QUSE9B5K`, attachment `SBW92NH4`. Its original verifier exits 0. The complete source note has 42 numbered equations, 275 unnumbered occurrences, and two figure crops.
- The two specialist reports under `numerical/ingestion/parallel/` and `numerical/ingestion/aliasing/` were read completely as judgment evidence; they are not substitutes for the canonical papers and are not wiki citations.

The source notes deliberately retain printed anomalies. This judgment does not convert the 2019 type-2 FFT-sign oddity, slash grouping, spherical-grid count issue, or the 2021 endpoint/saddle/sign/domain/logarithm issues into wiki formulas. The 2021 theorem must be printed as the source states it, followed by a caution that its stated `γ∈(0,1)` range is broader than the real-decay interpretation in the proof, where the denominator suggests `γ(2σ−1)>1`; this is a source/proof-regime caution, not a silent correction.

### Manual pages actually inspected

The source is the official v2.5.1 release tree at commit `679d9ae59fe0146c50da360e5a3fee70ae3aa646`, package 2.5.1. Its `conf.py` and built display still say 2.5.0; both facts must remain visible. The mirror contains all 35 content documents, but coverage is not the same as substantive reading.

Read in full:

- `FINUFFT_Overview`
- `FINUFFT_Transform_Definitions`
- `FINUFFT_CPU_Options`
- `FINUFFT_Status_Codes`
- `FINUFFT_Troubleshooting`
- `FINUFFT_GPU_Implementation`
- `FINUFFT_Inverse_Transform`
- `FINUFFT_Continuous_Transforms`
- `FINUFFT_Fourier_Evaluation`
- `FINUFFT_Real_Interpolation`
- `FINUFFT_Gaussian_Fields`
- `FINUFFT_Poisson_Python`

Read in bounded part:

- `FINUFFT_Python_Interface`: staged lines 170–230 and 848–1080, covering vectorized/Plan examples plus `Plan`, `setpts`, `execute`, and `execute_adjoint`. The middle simple multidimensional API entries were not read.
- `FINUFFT_C_Interface`: staged lines 468–705, the complete guru-plan interface. The simple/vectorized bodies were not read.
- `FINUFFT_Changelog`: staged lines 1–125, covering v2.5.1 and v2.5.0.

`FINUFFT_CPU_Options` was read completely but contributes no materially distinct selected claim. It is therefore recorded as an inspected raw-build-only exclusion, maps to no wiki page, and should not receive an `ingested` date from this judgment.

The remaining 20 content documents are **raw-build-only exclusions**, not inferred evidence: Acknowledgments, GPU C Interface, C Examples, GPU Migration, Developer Notes, Directories, Fortran Interface, Installation, GPU Installation, Julia Interface, MATLAB Interface, GPU MATLAB Interface, NFFT3 Migration, Performance, GPU Python Interface, References, Related Packages, Tutorials, Poisson MATLAB, and Users. The complete machine-readable source-page map records every one. The manual root README is navigation/provenance; it supports the FINUFFT entity's snapshot identity but is not automatically “ingested,” and no mirrored page receives an `ingested` date merely for existing.

## Duplicate decisions

The required two passes are complete.

First, actual `python3 99_SYSTEM/scripts/wiki_check.py find "<name>"` output exists for every adopted exact page name. The tool is alias-aware. All five proposed new names are unresolved; all seven existing names resolve to the exact page paths above. Supplemental alias calls also completed. Notable results:

- `FFT`, `Discrete Fourier transform`, `Fourier matrix`, and `Cooley-Tukey` resolve to `Fast_Fourier_Transform`.
- `Numerical integration`, `Trapezoidal rule`, and `Simpson's rule` resolve to `Numerical_Quadrature`.
- `Least squares` and `Best approximation` are ambiguous between the existing wiki page and a source-layer lecture note, while the exact basename `Least_Squares` is unambiguous. Preserve the existing aliases; new links should use the canonical basename.
- The proposed aliases for the five new pages were unresolved in the completed supplemental calls.

Second, I read all six candidate category READMEs. They establish these semantic boundaries:

- `12_LinearAlgebra_README`: no NUFFT, ES-kernel, NUFFT-aliasing, or parallel-spreading owner exists. Generic adjoints already belong to `Inner_Product_Space`; do not create `Adjoint_Operator`.
- `00_Entities_README`: no FINUFFT owner exists; create one flat `type: method` entity with no `up` field.
- `04_Calculus_README`: `Numerical_Quadrature` owns the Fourier-integral bridge; do not create `Fourier_Quadrature`.
- `12_Waves_and_Acoustics_README`: `Harmonics_and_Fourier_Series` owns off-grid Fourier-series evaluation; do not create `Fourier_Interpolation`.
- `05_Electrostatics_README`: `Poisson_Laplace_Equations` owns the periodic spectral-solve example.
- `10_Models_README`: `Gaussian_Random_Field_Media` can take a clearly separated stationary spectral-sampling implementation note, but remains primarily a level-cut-media page.

`duplicate-find-receipts.json` retains the commands, outputs, evidence files/session, alias results, and semantic decisions.

Rejected splits are deliberate: no `NUFFT_Conditioning` page (the source-backed material is bounded to NUFFT error caveats and the inverse problem on `Least_Squares`), no `Adjoint_Operator`, no `Fourier_Quadrature`, no `Fourier_Interpolation`, no `Toeplitz_Matrix` or `Conjugate_Gradient` from one tutorial, and no page per language/API/manual heading. Uninspected API sections and navigation/reference lists remain raw-build coverage only.

## Page specifications

### 1. New `Nonuniform_Fourier_Transform`

- Metadata: `type: concept`, `status: 🌿 growing`, tags `concept, linear_algebra`; aliases `NUFFT`, `Nonuniform FFT`, `Nonuniform fast Fourier transform`, `Non-equispaced fast Fourier transform`.
- Relations: `up: [[Linear_Algebra]]`; prerequisites `[[Fast_Fourier_Transform]]`, `[[Inner_Product_Space]]`; related `[[Exponential_Semicircle_Kernel]]`, `[[NUFFT_Aliasing_Error]]`, `[[Parallel_Spreading]]`, `[[FINUFFT]]`, `[[Numerical_Quadrature]]`.
- Index summary: “Three NUFFT types, sign/mode conventions, adjoints versus inverses, and spread-FFT-deconvolve complexity.”
- Claims: define types 1/2/3 with an explicit, internally consistent sign convention and centered even/odd mode sets; distinguish type-1/2 periodic coordinates from arbitrary type-3 sources/targets; explain direct `O(MN)` versus spread–FFT–deconvolve; identify the space-frequency extent controlling type-3 cost; and state that an adjoint is not an inverse.
- Exact sources: Barnett et al. 2019 §§1 and 3.1–3.3; `FINUFFT_Transform_Definitions` whole substantive page; current `FINUFFT_C_Interface` guru `setpts` and `execute_adjoint`; `FINUFFT_Python_Interface#^src-693fb6214ce497ea`; `FINUFFT_Overview#^src-b61f1ce2d53c4f18` and `#^src-bb5afcf91acb14fa`.
- Critical qualification: pinned C guru docs accept any real type-1/2 coordinate and fold it periodically, with rounding risk for huge magnitudes. This current API statement overrides the stale narrower coordinate-range prose in the math page.
- Exclude inverse algorithms beyond a link to `Least_Squares`, universal relative-error claims, current CPU kernel guarantees, performance marketing, and the 2019 printed type-2 sign anomaly.

### 2. New `Exponential_Semicircle_Kernel`

- Metadata: `type: concept`, `status: 🌿 growing`; aliases `Exponential of semicircle kernel`, `ES kernel`.
- Relations: `up` and prerequisite `[[Nonuniform_Fourier_Transform]]`; related `[[NUFFT_Aliasing_Error]]`, `[[Parallel_Spreading]]`, `[[FINUFFT]]`.
- Index summary: “Compact exponential-of-semicircle spreading kernel, its efficient evaluation, and the qualified 2019/2021 error story.”
- Claims: give the normalized compact-support kernel and separately its grid-width rescaling; distinguish two 2019 computations—numerical quadrature for the kernel Fourier transform/deconvolution factors versus generated piecewise-polynomial/Horner evaluation of spatial-kernel values in spreading/interpolation loops; distinguish the 2019 error statement from the 2021 proof; qualify ES/KB/PSWF asymptotic relationships; and warn that v2.5.1 CPU defaults to PSWF while GPU remains ES.
- Exact sources: Barnett et al. 2019 §§1.2, 3.1.1, 4.1–4.2, and 5; Barnett 2021 §§1, 3–5; `FINUFFT_Changelog#^src-13b56fe21738f5c7`; `FINUFFT_Overview#^src-68079183c36c65e1`.
- Exclude a finite-width universal constant, repaired paper formulas, or transfer of ES theory to current CPU PSWF.

### 3. New `NUFFT_Aliasing_Error`

- Metadata: `type: concept`, `status: 🌿 growing`; aliases `NUFFT aliasing error`, `NUFFT error bound`.
- Relations: `up: [[Nonuniform_Fourier_Transform]]`; prerequisites `[[Nonuniform_Fourier_Transform]]`, `[[Exponential_Semicircle_Kernel]]`; related `[[FINUFFT]]`, `[[Parallel_Spreading]]`.
- Index summary: “Absolute aliasing bounds for gridding NUFFTs, including the ES theorem's exact hypotheses and limits.”
- Claims: present the general error-matrix/tail representation and absolute `ℓ1→ℓ∞` bound; state the printed 2021 theorem only under its exact one-dimensional, type-1/2, exact-arithmetic, fixed-parameter, `w→∞` hypotheses; explain endpoint/saddle and phased-sinc cancellation; and list what the theorem does not establish.
- Exact sources: Barnett et al. 2019 §4.1; Barnett 2021 §§1–4, especially Theorem 1, general aliasing error, Fourier-transform asymptotics, and phased-sinc proof; current kernel boundary from `FINUFFT_Changelog#^src-13b56fe21738f5c7`.
- Exclude universal relative error, finite-width certification, multidimensional/floating/type-3 extrapolation, CPU-PSWF transfer, and silent anomaly repair.

### 4. New `Parallel_Spreading`

- Metadata: `type: concept`, `status: 🌿 growing`; aliases `Parallel NUFFT spreading`, `Parallel gridding`.
- Relations: `up` and prerequisite `[[Nonuniform_Fourier_Transform]]`; related `[[Exponential_Semicircle_Kernel]]`, `[[FINUFFT]]`.
- Index summary: “CPU subproblem and GPU output-driven strategies for race-aware, cache-aware NUFFT spreading.”
- Claims: bin sort and cache locality; CPU local subproblems, dynamic scheduling for clusters, race avoidance, and controlled global accumulation; GPU output-tile blocks, shared-memory staging, inner parallelism, avoided local atomics, and amortized global writes; architecture-specific implementations of one spreading atom.
- Exact sources: Barnett et al. 2019 §5; `FINUFFT_GPU_Implementation#^src-341a84e3cbc08314`, `#^src-374cd2e5893b45b8`, and `#^src-a66c3c22a552b289`.
- The GPU source is a living v2.4.1-era implementation description carried in the accepted v2.5.1 manual; describe it as implementation documentation subject to change, not a stable ABI or universal architecture rule.
- Exclude benchmark tables, hardware-independent speedups, ABI promises, and malformed upstream math roles.

### 5. New entity `FINUFFT`

- Metadata: `type: method`, `status: 🌿 growing`, tags `entity, software, finufft`; aliases `Flatiron Institute NUFFT`, `cuFINUFFT`; omit `up` and prerequisite fields under entity conventions.
- Related: the four new NUFFT concept pages plus `[[Fast_Fourier_Transform]]`.
- Index summary: “Pinned FINUFFT v2.5.1 software entity: transform scope, reusable plans, kernel split, tolerance semantics, and operating caveats.”
- Claims: official tag/commit/package and stale display-release distinction; supported types/dimensions and interface families; Plan/makeplan–setpts–execute/execute_adjoint–destroy lifecycle and reuse; C pointer-retention lifetime; moving points require `setpts`; CPU PSWF/GPU ES split; status/tolerance semantics; roundoff/conditioning/repeatability/numerical-adjoint cautions; setup-dominated small problems and possible direct/BLAS alternative.
- Exact sources: canonical manual README provenance; `FINUFFT_Overview` substantive blocks; `FINUFFT_C_Interface#^src-e20db61f1af66a8c`; `FINUFFT_Python_Interface#^src-4d8097549f5361ee`, `#^src-0e958f0977953bcc`, `#^src-61cc5516c3961ef4`, `#^src-693fb6214ce497ea`; `FINUFFT_Status_Codes#^src-7937087c7ee61b15`; `FINUFFT_Troubleshooting#^src-c762203609afde2f` and `#^src-c43a947c2123045d`; `FINUFFT_Changelog#^src-13b56fe21738f5c7`; the 2019 paper only for historical ES-library context.
- `ier=0` is completion success, not an independent accuracy proof; `ier=1` warns that the request was below supported tolerance and a best-effort transform was still produced.
- Exclude copied API catalogs, universal speed claims, verified-RTD-route claims, CPU ES guarantees, or “status zero certifies error.”

### 6. Add to `Fast_Fourier_Transform`

Append one short, attributed bridge: gridding NUFFTs spread/interpolate between nonuniform points and an oversampled regular grid, call an ordinary FFT, then deconvolve. Cite Barnett et al. 2019 §§3.1–3.2. Add related links to `Nonuniform_Fourier_Transform` and `FINUFFT`; preserve the complete existing page, sources, aliases, and index entry. Do not add throughput claims.

For this page and all six existing pages below, preserve the current category-README line verbatim. `page-claim-source-map.json` records that exact one-line summary and the complete intended final `up`/prerequisite/related relationships.

### 7. Add to `Least_Squares`

Append a bounded “inverse NUFFT as least squares” section: `Af≈y`, type 2 applies `A`, matching type 1 applies `A*`, and iterative normal-equation solves avoid forming the dense matrix. The 1D tutorial exploits Toeplitz structure; sampling gaps can make the problem ill-conditioned, and normal equations square the condition number. Cite `FINUFFT_Inverse_Transform#^src-c7426c8c79176c85`, `#^src-303b027e9e65b3a1`, and `#^src-1f17503ee88febe6`. Preserve all existing material and source ownership. Do not claim a general preconditioner or regularization scheme, and do not create `Toeplitz_Matrix` or `Conjugate_Gradient` from this single tutorial.

### 8. Add to `Inner_Product_Space`

Append one concrete Fourier-sampling example: under complex Euclidean inner products, taking the adjoint reverses spaces and conjugates the exponential phase; it is not inversion without additional unitary/normalization conditions. Numerical adjoint checks allow floating error and need not be bitwise repeatable. Cite the C and Python `execute_adjoint` blocks, `FINUFFT_Troubleshooting#^src-c762203609afde2f`, and the inverse tutorial. Do not create `Adjoint_Operator`; this page remains the generic owner.

### 9. Add to `Numerical_Quadrature`

Append “Fourier integrals as weighted exponential sums”: arbitrary target frequencies yield type 3; a regular frequency grid can be rescaled to type 1; quadrature/truncation convergence is separate from NUFFT `eps`; and the NUFFT accepts, but does not choose or certify, problem-appropriate adaptive/singular quadrature nodes and weights. Cite `FINUFFT_Continuous_Transforms#^src-a60be94374ff31e3`, `#^src-00a625b0e7b4ee08`, `#^src-0f925fd4923a243c`, `#^src-e4970d5554bd5333`, and `#^src-821a232939d9cff0`. Do not reproduce the tutorial's inconsistent substitution prose (`2y f(y)` after `x=y²`, while code evaluates `f` at squared nodes), and do not create `Fourier_Quadrature`.

### 10. Add to `Harmonics_and_Fourier_Series`

Append off-grid evaluation: type 2 evaluates known centered Fourier coefficients at arbitrary points, with period `L` mapped by `2πx/L`. The real-valued shortcut is only the demonstrated odd-`N` construction; in higher dimensions the saving applies along one dimension, while the implementation remains complex internally. Cite `FINUFFT_Fourier_Evaluation#^src-054e04aaa047f31a`, `#^src-3373dbcc24c07928`, and `FINUFFT_Real_Interpolation#^src-3068f5654f800c35`, `#^src-73bd7c1e6f148a71`. Do not extrapolate to even `N`, copy timings, or create `Fourier_Interpolation`.

### 11. Add to `Poisson_Laplace_Equations`

Append one explicitly periodic numerical route: `-Δu=f` requires zero mean; the zero Fourier mode fixes the arbitrary constant; a **known weighted non-Cartesian quadrature** feeds a type-1 coefficient sum, nonzero modes are divided by `|k|²`, and type 2 evaluates the solution at targets. Fully local adaptivity loses some advantage because the global Fourier truncation follows the highest required frequency. Cite `FINUFFT_Poisson_Python#^src-2709ac3f895ee3f8`. This does not cover arbitrary unweighted samples, free-space/nonperiodic boundaries, or any application energy/gradient derivation.

### 12. Add to `Gaussian_Random_Field_Media`

Append a clearly separate implementation note limited to **stationary** spectral models: a regular frequency quadrature grid, independent Gaussian coefficients scaled according to the PSD and weights, and type-2 evaluation at arbitrary spatial targets. Keep truncation, frequency quadrature, and NUFFT tolerance as separate error budgets. Cite `FINUFFT_Gaussian_Fields#^src-3b8611e01f621e90`. Do not recast every GRF as stationary, attribute the page's existing level-cut construction to the manual, or copy the sample's `x/P` FINUFFT coordinate call; that conflicts with the explicit `2πx/L` scaling in the inspected Fourier-evaluation tutorial.

## Citation and ownership rules for every page worker

- Use the accepted paper notes for claims they establish. Do not imply independent reading of Dutt–Rokhlin, Fourmont, Slepian, Kaiser-Bessel, or other works merely cited by those papers; historical attribution must say it is reported by the inspected paper.
- Use the now-accepted canonical manual links and the checked block locators in the machine map. Canonical source metadata points to the corresponding pinned GitHub path at commit `679d9ae59fe0146c50da360e5a3fee70ae3aa646`; an unverified ReadTheDocs tag route is not source provenance.
- Existing pages are additive-only. Preserve every earlier source, alias, frontmatter field, body section, and index/list entry; append a small attributed section and the minimum source/related additions.
- Do not add an `ingested` date to unused manual pages. Only source pages supporting an actually accepted wiki claim qualify, and the conductor owns dates/metadata/log/index/queue.
- Keep every claim's qualifications next to the claim. Do not turn tutorial examples into measurements or universal crossovers.
- No application energy/gradient derivation, new raw-source acquisition, code experiment, or numerical experiment is authorized.

## Recommended dispatch order

1. Dispatch `Nonuniform_Fourier_Transform` and `FINUFFT` together; their file ownership is disjoint and all link targets are already known.
2. Dispatch `Exponential_Semicircle_Kernel` and `Parallel_Spreading` next.
3. Dispatch `NUFFT_Aliasing_Error`, then the seven additive pages in whatever grouping respects the coordinator's worker cap.
4. Verify the complete 12-page wave as one graph after every target exists.

The `depends_on_new_pages` arrays in the machine map express intended wikilink/prerequisite relationships and final-wave verification needs, **not hard scheduling dependencies**. Workers may create links to the known target paths before the linked page is written. The accepted canonical manual has removed the former audit wait.

## Files and machine checks

- `page-claim-source-map.json`: authoritative 12-page metadata, claim-level citations/locators, dependencies, exclusions, and complete paper/manual source-page routing, including every raw-build-only exclusion.
- `duplicate-find-receipts.json`: actual exact-name outputs, supplemental alias outcomes, and category README semantic decisions.

No canonical vault file, source build, mirror, backend, research source, experiment, index, log, or queue was modified by this judgment task.
