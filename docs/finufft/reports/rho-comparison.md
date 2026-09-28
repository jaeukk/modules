# Direct-sum versus FINUFFT collective coordinates

Fixed-coordinate comparison; no optimizations performed. Three seeds, three saved states each, five species, 80 constrained modes per species: 3,600 complex values per tolerance. Both backends evaluate the same loaded particle coordinates. FINUFFT use is checked (no fallback).

rho is the unnormalized sum of exp(+i k dot r) over 400 particles of each species. Error is the complex magnitude |rho_FINUFFT - rho_direct|. Four threads, original CPU affinity/environment. The original requested tolerance is 1e-12.

| Tolerance | States | Max absolute error | RMS absolute error | Max error / 400 |
|---|---|---:|---:|---:|
| 1e-12 | initial | 6.730632e-11 | 9.690884e-12 | 1.682658e-13 |
| 1e-12 | direct_ground | 6.979158e-11 | 1.026453e-11 | 1.744790e-13 |
| 1e-12 | finufft_ground | 4.972264e-11 | 9.551701e-12 | 1.243066e-13 |

An additional attempted 1e-14 tolerance test on the first initial state triggered the integration fallback; the harness rejected it explicitly. No tighter-tolerance accuracy claim is made. The partial TSV from that rejected test is preserved but excluded from statistics.

Ground-state constrained rho values are near zero: percentage errors relative to each rho are ill-conditioned and are not the primary metric. Normalizing by 400 is a particle-count normalization, not a relative error with respect to rho. The direct sum is a double-precision reference, not an exact-arithmetic oracle.

This isolates density evaluation from optimizer trajectory differences. Agreement of constrained collective coordinates does not imply identical particle configurations. Only the constrained modes at chi_i approximately 0.1 were tested here. Full modewise complex values, per-species statistics, source hashes, and build commands are retained beside this report. Production code and original simulation outputs were not modified; the sweep remains stopped.
