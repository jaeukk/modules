# Independent FINUFFT manual-mirror audit

**Verdict: DO NOT PUBLISH this snapshot.** The audit completed successfully and found repairable defects; that is not acceptance of the mirror. Eight of 35 documents fail a content-structure or offline-reference check; all 35 also need the coordinator-adopted source-URL provenance correction (F3). The other 27 pass the implemented content checks, subject to the shared publication prerequisites and stated limitations below.

Audit date: 2026-09-22. Auditor: Codex (GPT-6), independent of the converter. Owned writes: `manual-audit/` only. No source, converter, staged mirror, canonical vault, backend code, or numerical experiment was modified or executed. `snapshot/` freezes the 35 audited Markdown files; `results.json` records source-body and output hashes so this verdict cannot silently transfer to a repaired revision.

## Authority and exact coverage

Read the adopted SOFTWARE_MANUAL_ROUTE.md and MANUAL_DESTINATION.md; the runtime-resolved vault's actual 40_Manuals_README, full Meep mirror README, meep-doc-refresh skill, Wiki_Schema, 20_Notes_Policy, and `_shared/contracts/document-note.md`. The later software-mirror decision governs this audit over the older destination document's condensation wording. This is source-mirror fidelity testing, not a claim to have extracted or read 236 PDF pages. None of the 17 bounded PDF scopes is represented here as a completed extraction Task, and no PDF verifier/self-report pass was fabricated.

* The release archive's **own tar PAX comment** is `679d9ae59fe0146c50da360e5a3fee70ae3aa646`; this is stronger evidence than merely repeating a converter receipt. Its SHA-256 is `809aa60ea4bf11a9976642d30d6bd2634da0083bdebd580d64705bc4aecfcd0e`. All **619 regular archive members** match the extracted source tree byte for byte. The existing provenance ties that archive to official tag **v2.5.1**. This audit did not independently refetch the tag.
* Independent traversal of the **resolved HTML toctree** from index yields exactly the 35 manifest documents, once each: index/overview plus 27 other top-level documents and seven tutorials. There is no missing or duplicate content document, including GPU, performance, API/autodoc, or tutorial content.
* Traversing **all** body hyperlinks reaches 37 source documents because `overview` and `latexindex` are additional wrappers for already included material. Their content does not require two extra copied notes, but their references require alias routing: see F2. These are not silently discarded from the audit.
* Inspected real RST/include/autodoc inputs and resolved `.fjson` bodies, not just converter receipts. Source archive, `docs/conf.py`, module `__init__.py` files, relevant RST directives, and Sphinx log corroborate the resolved dependency route. `structure-results.json` records the direct include/literalinclude/autodoc dependencies for 12 documents.

## Independent content and structure results

The audit never imports or invokes `build_manual_mirror.py`. Beautiful Soup parses resolved Sphinx HTML; Pandoc independently parses the emitted Markdown as GFM with math and wikilinks. The use of Pandoc as a reader is a shared-tool limitation, not replay of the converter's HTML-to-Markdown algorithm. `--preserve-tabs` and `--wrap=none` prevent the audit renderer itself from introducing tab expansion or text wrapping.

| Checked source content | Independently derived inventory | Result |
|---|---:|---|
| Math | 339 total: 323 inline, 16 display; 7 printed labels | PASS: ordered parsed expressions, display/inline type, exact raw source TeX substrings, and source labels checked; no missing expressions |
| Code blocks | 264 | FAIL F1: 34 documents preserve the parsed code sequence byte for byte, modulo final newline removal; GPU Implementation parses into 5 blocks instead of 3, including 2 spurious prose blocks |
| Headings | 153 | PASS: ordered levels and full text agree in every document |
| Paragraph prose | 1,198 source paragraphs | Full ordered text comparison passes 34 documents; GPU Implementation fails when prose becomes code |
| List items | 637 source items | FAIL F1: item depth and complete item content differ in four documents; GPU Implementation loses two real nested list items to code |
| Definition associations | 97 term/description pairs across CPU/GPU Python APIs | PASS: full ordered term+description sequences, including signatures/defaults, remain contiguous; all 30 API IDs remain present and statically linkable |
| HTML tables | 0 | Verified directly across all 35 source bodies; no invented Markdown table and no unsupported table-fidelity claim |
| Images | 62 references to 62 unique files | PASS: ordered identity, alt text, and byte hashes agree with resolved files and original repository files; all are available offline |
| Source targets | 234 source IDs → 195 emitted block IDs | All 195 emitted anchors occur once as top-level blocks, with a preceding semantic block; 233 source-ID associations match the independent semantic-block comparison, one is uncertain (U1) |
| Text links | 474, including 279 mapped internal references | Every mapped target and anchor resolves statically; all source destinations agree after URL-encoding normalization, but 10 additional content references remain external (F2) |

Full content comparisons normalize prose whitespace, preserve character order and special symbols, and compare math and code separately. They are not count-only gates. Each source paragraph, definition term, and caption-like unit is also checked against the ordered output text; every list item is compared with its complete descendant content, which catches code/prose escaping from a list even when overall text and counts match. Full per-document inventories, differences, math bodies/labels, code counterexamples, and definition descriptions are in `results.json`; list trees and every source-ID association are in `structure-results.json`.

The seven HTML labels are index `(1)`, math `(1),(2),(3)`, continuous-transform tutorial `(1),(2)`, and inverse-transform tutorial `(1)`. These deliberately do not adopt the PDFs' chapter-based numbering. Raw TeX is present unchanged apart from documented delimiters and appended `\tag{}`; a whitespace difference introduced by the independent HTML writer inside math was resolved by checking the original Markdown, not reported as a source defect.

## Required repairs with concrete counterexamples

### F1 — Protected code restoration escapes list containers

The converter inserts protected code as multiline text with column-zero fences after Pandoc has already indented its placeholder for a list. Restoring the token does not preserve the list context. Counts alone miss most of this defect.

* **GPU Implementation** (`impl_gpu`), frozen `FINUFFT_GPU_Implementation__STAGED.md` lines 65–105: the Shared Memory code fence opens at column zero. At line 77, `The **Output-driven approach** does:` still has six leading spaces and therefore becomes an **indented code block** after the list has ended. Its two nested bullets appear literally as code; the paragraph beginning `The parallelism is flipped` becomes another code block. The actual source has three code blocks, not five, and 11 list items, not nine. See `sample-impl_gpu.png` and independently rendered supplement PDF pages 89–90 (printed 86–87).
* **GPU Migration** (`cufinufft_migration`), `FINUFFT_GPU_Migration__STAGED.md` lines 42–63: the first list item ends after `which previously had the signature`. Both signatures, `and now has the signature`, and the explanation of `int64_t`/`maxbatchsize` become sibling blocks outside that item. The same defect affects list items 2, 3, and 5. The source binds each before/after signature pair and explanation to its API-function item.
* **Developer Notes** (`devnotes`), `FINUFFT_Developer_Notes__STAGED.md` line 72 onward: the `Installing MWrap` list item loses its code block and the concluding `to your make.inc, and then you can use the make mex task.` from its container. Its text remains somewhere in the document, so full-text/count equality alone would wrongly pass it.
* **Troubleshooting** (`trouble`), `FINUFFT_Troubleshooting__STAGED.md` line 78 onward: CPU list item 9 and GPU list item 18 end before their commands/results. The CPU explanation beginning `Extract the relevant transform type` and GPU explanation beginning `Check if your time is dominated by the plan stage` no longer belong to the source bullet.

**Repair:** preserve each protected block's list nesting during serialization/restoration; then rerun full code sequence and list-subtree comparisons. Do not merely indent the first fence or suppress the extra parsed blocks in the audit.

### F2 — Ten source-owned internal cross-references require the network

All ten are listed in `structure-results.json.externalized_content_links`. Five occur in **Overview**, two in **Troubleshooting**, and one each in **GPU Installation**, **MATLAB Interface**, and **Continuous Transforms**.

Examples: Overview line 58 links equation `(1)` to `https://finufft.readthedocs.io/en/v2.5.1/overview/#equation-1d1`; GPU Installation line 36 sends `overview page` to `.../latexindex/#index`; MATLAB Interface line 54 and Continuous Transforms line 121 also send `(1)` to the online wrapper. The relevant equation, GPU heading, and “Do I even need a NUFFT?” heading already exist in the local Overview note. The source's duplicate overview-label warnings explain the wrapper destinations, but do not justify an offline mirror losing these local references.

**Repair:** define explicit `overview`/`latexindex` source aliases to the manifest Overview note and map their fragments to the equivalent local semantic blocks. Keep exactly 35 content notes. Recheck every source hyperlink destination and all block associations after changing routing.

### F3 — Replace unverified rendered-site provenance URLs with pinned official source URLs

The current metadata uses constructed ReadTheDocs routes such as `https://finufft.readthedocs.io/en/v2.5.1/index.html`. The archive and local Sphinx build are verified; an actually published RTD tag route was not verified, so these links must not imply that the mirror was acquired from that site. The coordinator explicitly adopted the conductor's correction during this audit: every `source_url` should identify its authoritative pinned GitHub source, for example `https://github.com/flatironinstitute/finufft/blob/679d9ae59fe0146c50da360e5a3fee70ae3aa646/docs/index.rst`, with the local resolved-build dependency provenance recorded separately. Use the corresponding real source path for every note. If the RTD URL is retained as a convenience field, label it as an unverified rendering route rather than the acquired source. Keep the initial snapshot frozen and apply this correction only in the repair revision. All 35 rows in the per-document map flag F3.

### U1 — One nearby target association is uncertain

Index source ID `id1` comes from the second `.. _index:` in included `overview.src` and precedes the three logos. The converter binds it to `src-81c3e56be2b75b70` after the long first prose paragraph, skipping the first actual semantic block (the logo group). This is not a naked anchor, but it is not the nearest faithful source block. Bind this source alias to the overview heading or logo block, or explicitly justify its intended front-page target. See the complete source/actual block text in `structure-results.json.summary.semantic_mismatches`.

The 39 other source-ID aliases coalesced onto existing semantic targets are adjacent Sphinx section/label aliases, not 39 duplicate emitted block-ID lines. No emitted block has multiple competing ID lines. Their mapping is recorded transparently; this audit treats equivalent source aliases as one semantic target rather than inventing duplicate rendered blocks.

## Metadata, source spelling, and publication prerequisites

Except for the source-URL provenance correction in F3, all 35 notes correctly record the tag, commit, primary `docs/*.rst` path, `verbatim: true`, converter attribution, and the honest **package/tag 2.5.1 versus stale display 2.5.0** distinction. No `ingested`, invented DOI, citekey, or Zotero identity is present. `agent: deterministic Sphinx/Pandoc mirror via Codex, 2026-09-22` truthfully attributes the method; it does not claim a PDF reader completed every page. Canonical note basenames contain no spaces/hyphens, have at most four underscore-separated words, and collide with no existing Markdown basename in the runtime-resolved vault.

**Dependency provenance is still pending for publication.** A primary RST path alone is incomplete. Concrete examples: Overview includes `docs/overview.src`; Python and GPU Python resolve package docstrings; C includes `c1d.doc`, `c2d.doc`, `c3d.doc`, `cguru.doc`; MATLAB includes `matlabhelp.doc` and a literal `Contents.m`; Changelog includes `CHANGELOG`. Record the direct dependency map and pinned build environment (Sphinx 9.1.0, documented autodoc mocks/source import paths, dependency lock and resolved body identity) in the published provenance/manifest. Twelve affected documents are marked PENDING in the per-document map. No claim is made that the primary source path itself is false.

**Pending navigation/assets owned by conductor:** the canonical README and parent navigation, complete mechanical rewrite ledger, official source/archive links, original 146-page PDF and separate 90-page supplement, repository LICENSE and NOTICE with attribution. These are absent from this staged content-only mirror and were explicitly assigned to conductor, so their absence is a publication prerequisite, not an invented content-extraction failure. Preserve the source's 2017–2026 Simons Foundation attribution. Canonical link/alias checks, final wiki ingestion, selective `ingested` dates, log/review tasks, and final publication remain conductor-owned.

Source-owned spelling is preserved: e.g. `emphasise` in mathematical definitions, `specifed` in API parameter prose, and `We in 2025 we added` in GPU MATLAB were not silently corrected. GPU Implementation's three malformed `:math:` roles already appear literally in resolved Sphinx HTML and the supplement PDF page 90; they are upstream defects, not lost source math nodes. Keep/disclose them instead of silently inventing corrected math. No unresolved converter placeholders, include/autodoc directives, unavailable images, or root-escaping relative content links were found. GFM auto-links some formerly bare URLs and internal API link labels lose their code font while retaining their text; these format changes should be enumerated in the final rewrite ledger.

## Actually viewed visual evidence

Built local HTML from the independently parsed Markdown and rendered it with headless Chromium plus local MathJax 3; all screenshots below were actually opened with the image viewer. Standalone block markers are hidden only for these visual samples; their original text and association are tested separately. This is a browser rendering approximation, **not an actual Obsidian runtime click test**.

| Viewed Markdown/HTML sample | Viewed independent PDF reference | Observation |
|---|---|---|
| `sample-math.png`, full mathematical-definitions page | `reference-original-22.png`, original PDF page 22 / printed 18 | Readable piecewise sets and type-1/type-2 formulas; HTML labels (1),(2) correctly differ from PDF (4.1),(4.2). Sample also renders type 3, but the corresponding next PDF page was not viewed. |
| `sample-python.png`, nufft1d1 API definition | `reference-original-106.png`, original page 106 / printed 102 | Signature/defaults, literal formula code, parameter association and source typo preserved. Flat bold terms are readable; source semantic nesting was tested for all 97 associations. |
| `sample-performance.png`, introduction and first four 1D type-1 figures | `reference-supplement-17.png`, supplement page 17 / printed 14 | First chart's three versions, stacked components, axes and 1.00x/1.03x labels agree; four local images load. No new performance computation was run. |
| `sample-impl_gpu.png`, full GPU implementation page | `reference-supplement-89.png` and `reference-supplement-90.png`, pages 89–90 / printed 86–87 | F1 is directly visible: nested prose and bullets become code in Markdown but remain prose/lists in PDF. The malformed math-role text is literal in both upstream PDF and Markdown. |

Only **five distinct PDF pages** were rendered/viewed, at 160 dpi for legible whole-page structure; none was illegible and no missing equation was reconstructed. The samples are bounded structure checks, not full-PDF reading. `visual-render-results.json` records 33 rendered math nodes without MathJax errors on the math sample and successful loading of the four performance images. Original/supplement page counts were checked with PyMuPDF. The other 231 PDF pages were not visually read.

## Limitations and rerun instructions

* Full-machine coverage is against the resolved Sphinx content, not every possible alternate rendering of source RST. Original archive equality and real input inspection establish source identity; they are not standalone fidelity acceptance.
* Prose comparison ignores whitespace; character-level content and order are checked, but arbitrary typography is not pixel-matched. Math raw substrings plus parsed expressions and exact code comparisons are stronger than the prose normalization. Source CSS layout, image percentage sizing, and admonition styling are mechanically reduced and require disclosure.
* Source-ID associations are independently compared, but receipt mappings are read as the proposed mapping under test. Hash-derived IDs are not regenerated from the converter. No Obsidian mobile/desktop renderer or actual anchor click was tested; syntax and nearest-block attachment are static checks.
* The visual sample is four Markdown structures and five PDF pages. Only the first performance figure is compared to a PDF page; all 62 assets nevertheless have independently matching bytes and ordered alt/association checks.
* Source links still pointing at upstream `master` are source-owned external references, preserved as written. The snapshot metadata pins the manual; this audit does not promise those external destinations stay unchanged or remain reachable offline.
* `audit.py`, `structure_checks.py`, and `link_checks.py` are evidence generators; a successful process exit means the evidence was generated, **not that publication passed**. Read their explicit findings and this verdict. They write only under this audit directory. `render_samples.py` and `render_samples.cjs` reproduce the bounded visuals; the latter contains local executable paths and is an audit tool, not a vault dependency.

Conductor should repair F1/F2/F3, settle U1, add publication provenance/navigation/license/PDF assets, and rerun independent body/structure/link checks plus the affected visual sample before canonical publication. Do not clear this verdict using only converter hashes or matching inventory counts.

## Per-document map

PASS means the implemented content checks passed, not final publication acceptance. PENDING dependency provenance is an explicit unfinished metadata requirement. U1 is the one uncertain nearby anchor association.

| Source document | Content | Findings | Math / code / images | Dependency provenance | Source URL | Anchor association |
|---|---|---|---|---|---|---|
| `ackn` | PASS | — | 0 / 0 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `c` | PASS | — | 0 / 15 / 0 | PENDING | REPAIR F3 | PASS (static) |
| `c_gpu` | PASS | — | 3 / 18 / 0 | PENDING | REPAIR F3 | PASS (static) |
| `cex` | PASS | — | 11 / 12 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `changelog` | PASS | — | 0 / 1 / 0 | PENDING | REPAIR F3 | PASS (static) |
| `cufinufft_migration` | FAIL | F1 | 0 / 8 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `devnotes` | FAIL | F1 | 0 / 2 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `dirs` | PASS | — | 0 / 0 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `error` | PASS | — | 2 / 1 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `fortran` | PASS | — | 5 / 6 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `impl_gpu` | FAIL | F1 | 0 / 3 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `index` | FAIL | F2 | 36 / 0 / 5 | PENDING | REPAIR F3 | UNCERTAIN U1 |
| `install` | PASS | — | 0 / 35 / 0 | PENDING | REPAIR F3 | PASS (static) |
| `install_gpu` | FAIL | F2 | 0 / 11 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `julia` | PASS | — | 0 / 0 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `math` | PASS | — | 33 / 0 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `matlab` | FAIL | F2 | 2 / 16 / 0 | PENDING | REPAIR F3 | PASS (static) |
| `matlab_gpu` | PASS | — | 2 / 12 / 0 | PENDING | REPAIR F3 | PASS (static) |
| `nfft_migr` | PASS | — | 5 / 4 / 0 | PENDING | REPAIR F3 | PASS (static) |
| `opts` | PASS | — | 13 / 4 / 0 | PENDING | REPAIR F3 | PASS (static) |
| `performance` | PASS | — | 0 / 0 / 38 | primary source correct | REPAIR F3 | PASS (static) |
| `python` | PASS | — | 0 / 33 / 0 | PENDING | REPAIR F3 | PASS (static) |
| `python_gpu` | PASS | — | 0 / 23 / 0 | PENDING | REPAIR F3 | PASS (static) |
| `refs` | PASS | — | 2 / 0 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `related` | PASS | — | 1 / 0 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `trouble` | FAIL | F1,F2 | 40 / 5 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `tut` | PASS | — | 0 / 0 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `tutorial/contft` | FAIL | F2 | 73 / 10 / 6 | primary source correct | REPAIR F3 | PASS (static) |
| `tutorial/grf` | PASS | — | 7 / 1 / 1 | PENDING | REPAIR F3 | PASS (static) |
| `tutorial/inv1d2` | PASS | — | 51 / 11 / 2 | primary source correct | REPAIR F3 | PASS (static) |
| `tutorial/peripois2d` | PASS | — | 21 / 8 / 4 | primary source correct | REPAIR F3 | PASS (static) |
| `tutorial/peripois2d_python` | PASS | — | 21 / 11 / 4 | primary source correct | REPAIR F3 | PASS (static) |
| `tutorial/realinterp1d` | PASS | — | 2 / 2 / 0 | primary source correct | REPAIR F3 | PASS (static) |
| `tutorial/serieseval` | PASS | — | 7 / 12 / 2 | primary source correct | REPAIR F3 | PASS (static) |
| `users` | PASS | — | 2 / 0 / 0 | primary source correct | REPAIR F3 | PASS (static) |
