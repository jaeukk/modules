# Independent FINUFFT repaired-mirror re-audit

**Verdict: ACCEPT the revision 2 content and fully prepared publication package, pending the conductor’s canonical copy/readback/link/alias gates.** No remaining content or provenance defect was found by this audit. This is acceptance of reviewed bytes, not a claim that canonical publication or wiki ingestion has occurred.

Audit date: September 22, 2026. Owned writes: `manual-reaudit/` only. No converter, staging, source, initial audit, canonical vault, backend, or numerical experiment was modified or executed. The publisher was exercised only in an owned synthetic vault. No subworkers were dispatched.

## Revision identity

The task began against manifest `b08fa143c5a55387a93f0c0b2f9895f3a63e3c58c70ec4b411f77969f301cf7e`. During review, coordinator message `msg_88a1e10270ab` authorized a narrow format amendment after the real vault checker misparsed three literal-bracket/wikilink adjacencies. I preserved that complete repair as `revision1/`, paused staging reads, notified the conductor, and resumed after coordinator message `msg_794c62354014` confirmed the new freeze. This verdict binds **only revision 2**, independently rerun across all 35 documents. `revision2/` preserves the final notes, receipts, assets and frozen package artifacts.

| Final frozen artifact | SHA-256 |
|---|---|
| `manual-mirror/manifest.json` | `34b8cbf104ae33436e6b7b71f633c3465d71f8739aa455cfb650312dd9e1e9e3` |
| `build_manual_mirror.py` | `d877fda5c709bbe759a860735a9522f2251381b6e6e13693268f1228d863e6df` |
| `Flatiron_FINUFFT_2026_README__DRAFT.md` | `1005046c2f17eef794b06e1ea4533fde7afeb960a999d1d0e53f0ea1d31e93ee` |
| `acceptance/manual-canonical-provenance.json` | `c117bf569cc681a6f736df83a79fee036f61ef63147a972ad334e6ab899e96f1` |
| `publish_manual.py` | `c014733eaa4f65231b4dac4133bed67d3cb9eb7344619a0b071f612d45ea38f1` |

`fullcheck-results.json` is the machine-readable verdict, complete per-document acceptance map, evidence-hash ledger and explicit remaining publication prerequisites. `package-results.json` records exact input hashes and all 104 proposed publication destinations/bytes/hashes. The initial audit’s rejection remains intact and its report hash is verified in the provenance sidecar.

## Authority and method

Read BRIEF, numerical/REPORT, SOFTWARE_MANUAL_ROUTE, MANUAL_DESTINATION, MIRROR_AUDIT_SPEC and the complete initial audit report. Resolved the actual vault from `~/.claude/zotero-obsidian.local` and read its Manuals README, full Meep mirror README, Policy, Wiki_Schema, meep-doc-refresh and wiki-ingest instructions, and document-note extraction contract as fidelity guidance. Policy/input hashes are recorded. The adopted software-manual route supersedes the older destination document’s PDF-condensation wording. This review does not fabricate PDF extraction passes or treat the 17 abandoned extraction scopes as completed.

The initial audit’s independent HTML/Markdown reader checks were copied into this owned directory and reviewed before rerunning; the source-URL expectation was updated and the link checker was made explicitly alias-aware. These checks parse actual resolved Sphinx HTML and independently parsed emitted Markdown; they never import or run the mirror converter. Additional newly written checks compare container ancestry and list types, inspect actual alias-fragment semantic targets, discover dependencies from pinned source files, validate the publication inputs, and compare the two repaired revisions. The conductor’s copied check exits were not used as acceptance evidence.

## Separate findings

- **F1 — ACCEPT repaired.** Every ordered code sequence agrees exactly modulo final newline removal. Full source list-item text/depth and separate ordered `ul/ol/li/pre/blockquote` ancestry agree across all 35 documents, including list start values. GPU Implementation has exactly three code blocks and 11 list items; prose and nested bullets no longer become code. Migration signatures/explanations stay within the correct API-function bullets; Developer Notes and Troubleshooting preserve their code/explanation containers.
- **F2 — ACCEPT repaired.** All 474 source text links match in order: 289 local references, including all ten former wrapper routes, and 185 source-owned external references. The ten `overview`/`latexindex` references are checked against their actual source fragment IDs and source target text, corresponding index semantic targets, and emitted local destinations. All local anchors exist exactly once. The old structure helper still calls these ten links “externalized” solely because its source-document name map excludes duplicate wrappers; its recorded warning is explicitly superseded by `link-results.json` and `package-results.json.alias_links`, not silently suppressed.
- **F3 — ACCEPT repaired.** All 35 source URLs point to the correct `docs/*.rst` file at pinned commit `679d9ae59fe0146c50da360e5a3fee70ae3aa646`. Each note describes the local Sphinx/Pandoc conversion, preserves `verbatim: true`, and distinguishes tag/package 2.5.1 from displayed 2.5.0. No invented DOI, citekey, Zotero identity or premature `ingested` date appears.
- **U1 — ACCEPT by semantic judgment.** Actual `overview.src` starts with the duplicate `.. _index:` and explicitly describes itself as front-page content without the header, immediately before the logos. The index RST supplies the whole-block H1 immediately above the include. Coalescing `id1` onto that H1 is faithful front-page navigation and no longer jumps past the logos to the long paragraph. The nearest-logo helper’s mismatch is preserved in evidence: 233 targets pass literal semantic-text comparison, while this one is explicitly accepted for its block meaning. All 234 IDs map to 194 unique emitted anchors; 30 API targets remain linkable.
- **Revision 2 bracket amendment — ACCEPT.** Exactly two CPU Python and one GPU Python return-type lines change `\[[[` to `&#91;[[`; all other note bytes equal revision 1. Full text, API associations, destinations, and source bracket display are preserved. The new rewrite is disclosed in the converter receipts, provenance and README. The permissive raw checker’s link-extraction pattern now finds no bracket-contaminated target.

## Complete coverage and source fidelity

| Independently checked inventory | Result |
|---|---|
| 35 navigation-reachable content documents | Exact one-to-one manifest coverage; duplicate wrappers route locally without extra notes |
| 339 math expressions (323 inline, 16 display), seven printed HTML labels | Ordered parsed forms and exact raw TeX bodies preserved |
| 264 code blocks, 153 headings, 1,198 paragraphs | Ordered content passes; prose whitespace normalization disclosed |
| 637 list items, 97 API definition associations | Full item content/depth and term/description continuity pass |
| 62 image references / 62 unique assets | Ordered identity/alt text and bytes match resolved and original source images |
| HTML tables | Zero in actual resolved bodies; no invented table-fidelity claim |
| Source archive | SHA-256 `809aa60ea4bf11a9976642d30d6bd2634da0083bdebd580d64705bc4aecfcd0e`; PAX commit matches; all 619 regular members match extracted sources |
| Note names/aliases | No planned basename or alias collision in 5,645 actual vault Markdown files; external metadata/README note targets resolve uniquely |

Every document has explicit ACCEPT, note hash and resolved-body hash in `fullcheck-results.json.per_document`. Code, mathematics, source misspellings and upstream malformed GPU `:math:` roles remain source-faithful. No placeholder, unresolved raw include/autodoc directive, missing image, or escaping relative body link was found. External upstream links remain source quotations of destinations, including `master`; they were not refetched or certified as version-matched evidence.

## Publication package actually reviewed

The prepared README lists all 35 notes, correct archive and both PDF hashes, 146/90 page counts and printed-folio offsets, source/installed versus stale display/module versions, selective-ingestion semantics, and the mechanical rewrite narrative. It does not claim all-PDF reading. Exact upstream LICENSE/NOTICE and Simons attribution are present in the publication plan.

The provenance’s 35 primary sources, resolved body hashes, note hashes and source-target maps match actual files. Independent line-by-line directive discovery agrees for all direct and recursive includes/literalincludes and autodoc inputs: 12 affected documents. Entire CPU/GPU Python package `.py` inventories, including files absent from the directly displayed API, are enumerated and hash-checked. Actual pinned `conf.py` corroborates import paths, mocks and stale display release; dependency lock, supplemental wrapper and portable configuration agree with build artifacts. Source archive/member equality independently backs the pinned-file identity; no upstream member was modified. This is verification of recorded build inputs and outputs, not a fresh Sphinx/LaTeX rebuild or an independently refetched release tag.

The plan contains **104 distinct outputs**: 35 notes, 62 images, README, archive, two PDFs, provenance JSON, LICENSE and NOTICE. All exist with reviewed bytes; no premature canonical copy was required. Synthetic-vault tests of the unchanged frozen publisher passed: conflicting last-planned NOTICE rejected before any copy, dry run wrote no vault file, all 104 files were exclusively created and read back, identical rerun passed, and unaccepted gate, mismatched manifest hash and changed note bytes were rejected.

Publisher limitation: assertions require normal Python, not `-O`; the publisher independently pins manifest/notes/images/README but does not independently pin provenance/PDF/archive/license inputs. Therefore conductor finalization must compare **every** publication receipt/readback hash with the accepted package inventory, not merely trust process exit. The intentional provenance append of this final report/hash must be recorded separately while preserving all substantive source/build records and the initial rejection. This is a finalization requirement, not missing prepared content.

Remaining conductor work: append accepted re-audit identity, record acceptance, copy/read back the reviewed package, link parent navigation, run actual canonical raw wikilink and alias checks, and perform the separately owned wiki selection/ingestion/date/log/review/queue steps. Acceptance here does not preempt those gates.

## Actually viewed visual scope

Viewed with the image tool: repaired GPU Implementation and GPU Migration whole-page browser renders, mathematical definitions, the `nufft1d1` API sample, the performance introduction/first four figures, and the final bracket-context comparison. Chromium rendered locally parsed GFM with local MathJax; this is an approximation, not an Obsidian runtime click test. GPU Implementation and Migration were compared with actual pinned RST/resolved content and supplement PDF pages **89–90** and **87–88**, respectively. Their nested prose, bullets, signatures and code associations are visibly repaired; source malformed GPU math roles remain literal in both renderings.

Also viewed original PDF pages **22** (math) and **106** (API), and supplement page **17** (first performance figure). HTML equation labels deliberately differ from chapter-numbered PDF labels. All 33 math nodes in the definitions sample rendered without MathJax errors; four local performance images loaded. The final bracket sample shows all three amended return-type lines and six nearby parameter contexts beside their source equivalents; literal brackets, link labels and surrounding text agree.

Exactly **seven distinct PDF pages were viewed**; the other 229 were not visually read. PDF references were rendered at 160 dpi except supplement page 90 at 140 dpi; these were legible structural samples, not equation recovery. Some additional generated images, including Troubleshooting and automatically located supplemental references, were not opened and are not claimed as viewed. The five earlier content samples remain byte-unchanged between the two repaired revisions, as checked independently; the bracket sample was built and viewed after revision 2.

## Evidence and limits

`results.json`, `structure-results.json`, `link-results.json`, `package-results.json`, `final-check-results.json`, `publisher-test-results.json`, visual results/images and the preserved revisions are durable evidence. `fullcheck-results.json` hashes the scripts/results and every claimed viewed image. The scripts’ successful exit is not the gate: `assemble_report.py` asserts the full substantive results, retains the old-helper limitations and explicit U1 judgment, then emits this accepting verdict.

Prose normalization ignores whitespace, shared Pandoc parsing may share parser limitations, and CSS/admonition/definition-list typography is mechanically reduced as disclosed. Static anchor associations and bounded browser views do not prove every Obsidian renderer’s behavior. No numerical experiment, tutorial execution, full-PDF extraction, live external-link check or canonical publication is claimed.
