# Release validation — September 28, 2026

Fresh build and `potential/tests/finufft_repeatability/run_repeatability.py`: PASS, 101 required checks, zero failures. Twelve optimizations: three seeds, two backends, two separate-process repeats. All endpoints pass independent physical criteria; every same-backend matched pair has identical endpoint files and optimization pass records. Macro-enabled and macro-disabled API tests and the 100-repeat fixed-coordinate probe pass. FINUFFT use and its one-thread production default are observed without overriding options.

Conditions: five components, N=2000, chi per component 80/798, total packing fraction 0.2, four outer OpenMP threads, one internal FINUFFT thread. Existing three-pass `RelaxStructure_NLOPT()` is unchanged.

| Seed | Vanilla seconds r0/r1 | FINUFFT seconds r0/r1 | Cross-backend mean displacement |
| --- | --- | --- | ---: |
| 20260928 | 6.729 / 6.970 | 6.664 / 6.317 | 0.4309 |
| 20260929 | 6.894 / 6.859 | 6.206 / 5.997 | 0.4107 |
| 20260930 | 6.059 / 6.312 | 5.472 / 5.433 | 0.2864 |

Median optimization time: vanilla 6.794 s; FINUFFT 6.102 s. These are measured optimization times with potentially different evaluation counts. Cross-backend positional agreement is explicitly not an acceptance claim.

Full local receipts, source/library/binary manifests, API output, and endpoints: `_multiagent/tasks/finufft-release-validation-20260928/` in the surrounding research workspace. The original report is included as [release regression snapshot](reports/release-regression.md).

The factor-only diagnostic entry points were retired, with recoverable source and binaries outside the code repository. `git diff --exit-code -- potential/StructureOptimization.cpp potential/StructureOptimization.h` confirms no production optimizer modifications. No broad sweep was restarted.

## Type 1 unit suite

`bash potential/tests/finufft_type1/run_tests.sh`: all variants passed.

- Direct-only build: 40960/40960 checks; no libfinufft dependency.
- FINUFFT-enabled build: 40982/40982 checks.
- FINUFFT plus `HAVETWOFIXEDPARTICLES`: 40982/40982 checks.

The suite covers density/energy/force agreement, eligible and fallback modes, cache behavior, dimensions and cell geometries, finite differences, and fixed-seed smoke checks. Smoke timings are informational, not release benchmark claims.
