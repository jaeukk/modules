**Review: cco-finufft-blocks-r1, 854eea9..930d18a**

**No blocking findings.**

Reviewed the supplied diff, cache lifecycle, block scheduling, fallback paths, contact re-binning, and GetCCO initialization flow. The GetCCO wrapper applies after random/input selection and before ground/MD branching; its only explicit configuration mutation is `SetCellSize`.

**FOLLOW-UP — stale open item**

- **Location:** `docs/finufft/perf/REPORT.md:353`
- **Trigger:** Reading open item 13 after commit 930d18a.
- **Impact:** It incorrectly states that GetCCO still uses coarse contact cells. The wrapper at `apps/CCOptimization/main.cpp:235` addresses this.
- **Suggested correction:** Remove or mark item 13 resolved.

**Commands and results**

| Command/check | Result |
|---|---|
| `make -C apps/CCOptimization USE_FINUFFT=1 -j6 test-finufft test-finufft-forces test-species-concurrent test-indexed-view test-fused-contact test-finufft-blocks` | Exit 0. All suites passed, including 153 block checks, 168 plan-cache checks, 2,337 species-concurrent checks, and 916 fused-contact checks. |
| `TMPDIR="$PWD/review-tmp" bash potential/tests/finufft_type1/run_tests.sh` | Exit 0. Direct: 40,977 checks; FINUFFT and fixed-particle variants: 40,999 each. `ALL VARIANTS PASSED`. |
| `make -C apps/CCOptimization -j6 all && python3 apps/CCOptimization/tests/test_cli.py` | Exit 0. CLI tests passed, including GetCCO byte-for-byte repeatability. |
| `OMP_THREAD_LIMIT=1 timeout 90s apps/CCOptimization/test_finufft_blocks.out` | All 153 checks passed; no deadlock. |
| `OMP_THREAD_LIMIT=2 OMP_DYNAMIC=true timeout 90s apps/CCOptimization/test_finufft_blocks.out` | All 153 checks passed; no deadlock. |
| Reviewer-only nested-region probe | Wrapped the mixed-species fused evaluation in an outer two-thread region after plans were prepared. All 153 checks passed with `OMP_MAX_ACTIVE_LEVELS=1`, and with `OMP_MAX_ACTIVE_LEVELS=4 OMP_THREAD_LIMIT=4`. |
| `make -C apps/CCOptimization USE_FINUFFT=1 bench_finufft.out`, then `./bench_finufft.out opt 3 2000 2 0.4 8 1 finufft` twice | Both runs completed 324 evaluations with hash `8783c1a687821511`. Identical final energy: `3.5888547570346691e-20`; Direct reevaluation: `3.5966457967178987e-20`. |

The passing checks cover block-order rho reduction, Type 2 block-count identity, cold AllForce, fused/separate equality, copy/assignment ownership, and injected makeplan/setpts failures. Reduced-team and nested probes found no scheduling failure. These runs are evidence on this build and runtime, not an exhaustive race-detector or supported-platform verification.