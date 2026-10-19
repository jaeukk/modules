# Opus review: disposition and fixes

Author: Codex. Reviewer: Claude Opus 5.5, dispatched read-only through an isolated
policy copy; the installed policy was not edited. Audit budget: two rounds.

## Round 1

The reviewer found no blocking issue in the C++ evaluation logic and reported
three blockers outside it. See `audit-opus-round1.md`.

- **B1 fixed:** the direct-only CLI test now matches the actual `backend not
  available` error. The CLI also prints its selected schedule from the getter;
  the test checks selecting B and resetting to A. Direct-only and FINUFFT checks
  pass. This was also identified in self-review while the audit was running.
- **B2 fixed:** a local `.gitignore` excludes `__pycache__/`. Bytecode will not be
  staged. Review documents are intentionally committed; raw benchmark outputs
  are written outside the repository.
- **B3 fixed:** the new API and helper have Doxygen parameter/return documentation;
  Energy, AllForce and fused-evaluation docs distinguish the new schedule and its
  repeatability limit. The benchmark runner functions have NumPy-style docstrings.
- **C.21 note fixed:** both scope guards explicitly prohibit copying and moving.
- **Benchmark harness self-review fix:** initialize `ProgramStart` and `TimeLimit`
  before using the shared optimizer, as the existing benchmark does. A tiny
  bounded NLopt validation run emits valid JSON and reports four internal threads
  for every species. This is a harness smoke check, not the A/B performance run.

The remaining non-blocking notes are not landing criteria violations:

- Derived potential callbacks can have their own parallel regions. GetMultiCCO
  creates exact base species; Option B preserves virtual behavior for derived
  classes rather than rewriting their algorithms. Its no-enclosing-team guarantee
  refers to the base species path. Derived callbacks remain responsible for their
  own nested parallelism.
- The focused fixtures do not force a multi-block density partition. The existing
  block implementation is unchanged; Option B executes every returned block in
  order. The performance cases will record actual backends and internal threads.
- Async thread creation failure propagates safely as an exception; it is not
  silently treated as a completed contact calculation.

## Executed validation before commit

- Focused C++ tests: `OMP_THREAD_LIMIT=1`, `2`, and `4`, passive waiting.
- AddressSanitizer on the new tests and `MultiStealthyPotential.cpp`; compilation
  uses `-Wall -Werror`. No sanitizer errors or compiler warnings.
- Direct-only and FINUFFT CLI builds; new schedule/reset/missing-backend tests and
  the existing `tests/test_cli.py`.
- Existing text/MSPB export, MSPB archive and CLI, and offline PairStat checks.
- A small bounded NLopt harness validation (Option B only, ten particles/species,
  one evaluation cap per phase). Actual size comparisons start after commit.

Round 2 reviews only these fixes and verifies that the blockers are resolved.

## Round 2 and final checks

Claude Opus round 2 reports no remaining blocking findings; B1, B2 and B3 are
resolved. Its report remains unchanged in `audit-opus-round2.md`.

The final AddressSanitizer rebuild of both changed C++ sources, using
`-Wall -Werror`, completed after the mechanical scope-guard changes and passed.
The final direct-only and FINUFFT CLI builds and both CLI suites also passed.
Compiler logs contain no warning or error diagnostics. The final focused tests
passed with effective thread limits one and four (the earlier two-thread run
also passed). The two permitted Opus audit cycles are complete.
