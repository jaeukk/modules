# Release validation

Recorded 2026-10-19T06:44:40.022969+00:00. Implementation base: 856d6d1.
The release changes the shared NLopt history default to 40 and retains Option A.
All new production settings and the scheduling benchmark use contact strength 100.
The archived strength-1 test remains explicitly labeled historical evidence.

Passed on the local AMD Ryzen 5 7500F / WSL2 host:

- Warning-free direct and FINUFFT builds with the repository's -Wall flags.
- `make USE_FINUFFT=1 test-sequential-species`, including thread limits 4 and 1:
  direct energy/forces, finite differences, batching, motion, caches, toggles,
  fallbacks, exceptions, and joining the contact worker.
- `tests/threading/test_cli.py` in both builds: A/B/reset selection, default 40,
  explicit 40 producing identical bounded endpoints, explicit 100, rejection of
  negative/zero/fractional/malformed/overflow history, and unavailable FINUFFT.
- `tests/test_cli.py` in both builds: export, resume, PairStat, GetCCO, repeatability.
- `make USE_FINUFFT=1 test-mspb-pack`: archive API, exact snapshot restoration,
  offline PairStat, optimizer export/resume, mismatch rejection.
- `tests/run_export_test.sh` and `tests/run_pairstat_test.sh`.
- Updated benchmark smoke runs: default 40 and explicit 100 recorded correctly,
  contact-strength metadata 100, direct checks passing.
- Portable archived history harness rebuilt and completed a small strength-100
  smoke run, with objective counts matching all three NLopt phases.
- HTML report: all three interactive chart metrics, complete JSON download,
  390-pixel layout without page overflow, no JavaScript errors.
- `git diff --check`.

Della's GCC 11 build additionally diagnosed the preexisting MINOP candidate-energy
declaration. The follow-up initializes it within each iteration. The CLI suite
now exercises zero-iteration identity and decreasing energy after five MINOP
iterations. Both direct and FINUFFT builds pass those checks. This does not change
the NLopt workflow used by the scheduling and history comparisons.

The two previous Opus audit reports are adjacent to this evidence directory.
They reviewed the scheduling implementation before 856d6d1, not the later default
and report changes or a new Della production protocol. No additional audit cycles
were used for this release.
