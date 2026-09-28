# Five-component FINUFFT study

> Audit correction: the sweep is stopped; the status below is the last pre-stop snapshot. Paired initial files are byte-identical, but four-thread same-backend repeats also produce large endpoint displacements. These displacements cannot be attributed specifically to FINUFFT. See [pairing and reproducibility audit](../PAIR_AUDIT.md).

Updated: 2026-09-27T06:16:58.317204+00:00

Status: **running**; phase: step3.

N=2000 (400/species), 2D, equal radii, number density 1, total packing fraction 0.2. Four OpenMP threads on distinct physical cores; passive waiting. u0=e0=1. Complete reciprocal shells; chi is per component, M/(2*399).

Optimization uses the existing three-pass L-BFGS routine, at most 2000 evaluations/pass and 600 s/configuration. Reported seconds are actual optimization wall times, not physical simulation time. Each endpoint is independently recomputed using direct density/force sums and all-pair contact checks.

Acceptance: total reference energy<=1e-11, max component S<=1e-10, max force component<=1e-7, max relative overlap<=1e-6. A passing endpoint and an optimizer success code are distinct recorded facts.

## Step 2: paired validation timings

| Seed | Vanilla seconds | FINUFFT seconds | Vanilla accepted | FINUFFT accepted |
|---|---:|---:|---|---|
| 20260928 | 6.459 | 9.864 | True | True |
| 20260929 | 6.608 | 9.447 | True | True |
| 20260930 | 5.772 | 8.772 | True | True |

Matched-particle differences after periodic wrapping and global translation:

| Seed | Mean displacement | RMS displacement | Maximum displacement |
|---|---:|---:|---:|
| 20260928 | 4.1632e-01 | 5.0289e-01 | 2.0569e+00 |
| 20260929 | 3.9514e-01 | 4.7248e-01 | 2.1836e+00 |
| 20260930 | 2.8683e-01 | 3.4376e-01 | 1.9469e+00 |

## Step 3: FINUFFT-only random-start sweep

Vanilla optimization: **not run** in Step 3. Direct summation is used for endpoint validation only.

| Target chi_i | Actual chi_i | Completed / 10 | Accepted | FINUFFT median seconds | FINUFFT min–max seconds |
|---|---:|---:|---:|---:|---|
| 0.10 | 0.100251 | 10 | 10 | 9.290 | 8.894–9.865 |
| 0.20 | 0.197995 | 10 | 10 | 12.227 | 11.455–13.439 |
| 0.30 | 0.300752 | 10 | 10 | 19.380 | 18.136–20.737 |
| 0.40 | 0.401003 | 3 | 3 | 45.927 | 44.862–46.505 |
| 0.45 | pending | 0 | — | — | — |
| 0.50 | pending | 0 | — | — | — |
| 0.55 | pending | 0 | — | — | — |

## Individual simulations

| Stage | Target chi | Seed | Backend | Seconds | CC energy | Contact energy | max S | max overlap | max force | Accepted | Pass codes |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---|---|
| 2 | 0.10 | 20260928 | direct | 6.459 | 2.077e-20 | 0.000e+00 | 4.936e-20 | 0.000e+00 | 1.346e-12 | True | 1,1,1 |
| 2 | 0.10 | 20260928 | finufft | 9.864 | 5.916e-20 | 0.000e+00 | 2.777e-19 | 0.000e+00 | 1.778e-12 | True | 1,1,1 |
| 2 | 0.10 | 20260929 | direct | 6.608 | 4.519e-20 | 0.000e+00 | 3.306e-19 | 0.000e+00 | 1.596e-12 | True | 1,1,1 |
| 2 | 0.10 | 20260929 | finufft | 9.447 | 1.196e-19 | 0.000e+00 | 2.753e-19 | 0.000e+00 | 2.146e-12 | True | 1,1,1 |
| 2 | 0.10 | 20260930 | direct | 5.772 | 5.729e-20 | 0.000e+00 | 1.482e-19 | 0.000e+00 | 1.405e-12 | True | 1,1,1 |
| 2 | 0.10 | 20260930 | finufft | 8.772 | 3.017e-20 | 0.000e+00 | 1.005e-19 | 0.000e+00 | 1.622e-12 | True | 1,1,1 |
| 3 | 0.10 | 20261001 | finufft | 9.271 | 6.826e-20 | 0.000e+00 | 3.054e-19 | 0.000e+00 | 2.094e-12 | True | 1,1,1 |
| 3 | 0.10 | 20261002 | finufft | 8.894 | 4.529e-20 | 0.000e+00 | 1.387e-19 | 0.000e+00 | 1.651e-12 | True | 1,1,1 |
| 3 | 0.10 | 20261003 | finufft | 9.309 | 6.188e-20 | 0.000e+00 | 2.362e-19 | 0.000e+00 | 2.170e-12 | True | 1,1,1 |
| 3 | 0.10 | 20261004 | finufft | 9.437 | 1.290e-20 | 0.000e+00 | 2.486e-20 | 0.000e+00 | 8.140e-13 | True | 1,1,1 |
| 3 | 0.10 | 20261005 | finufft | 9.192 | 5.738e-20 | 0.000e+00 | 1.686e-19 | 0.000e+00 | 1.675e-12 | True | 1,1,1 |
| 3 | 0.10 | 20261006 | finufft | 8.971 | 8.855e-20 | 0.000e+00 | 3.788e-19 | 0.000e+00 | 2.172e-12 | True | 1,1,1 |
| 3 | 0.10 | 20261007 | finufft | 9.499 | 2.765e-20 | 0.000e+00 | 1.691e-19 | 0.000e+00 | 1.893e-12 | True | 1,1,1 |
| 3 | 0.10 | 20261008 | finufft | 9.865 | 8.613e-21 | 0.000e+00 | 3.586e-20 | 0.000e+00 | 7.570e-13 | True | 1,1,1 |
| 3 | 0.10 | 20261009 | finufft | 9.619 | 5.967e-20 | 0.000e+00 | 1.444e-19 | 0.000e+00 | 2.183e-12 | True | 1,1,1 |
| 3 | 0.10 | 20261010 | finufft | 9.081 | 2.586e-20 | 0.000e+00 | 1.001e-19 | 0.000e+00 | 1.906e-12 | True | 1,1,1 |
| 3 | 0.20 | 20261001 | finufft | 12.315 | 2.765e-20 | 0.000e+00 | 7.700e-20 | 0.000e+00 | 1.964e-12 | True | 1,1,1 |
| 3 | 0.20 | 20261002 | finufft | 12.909 | 1.615e-20 | 0.000e+00 | 2.592e-20 | 0.000e+00 | 1.215e-12 | True | 1,1,1 |
| 3 | 0.20 | 20261003 | finufft | 12.262 | 3.796e-20 | 0.000e+00 | 9.312e-20 | 0.000e+00 | 1.640e-12 | True | 1,1,1 |
| 3 | 0.20 | 20261004 | finufft | 13.126 | 2.772e-20 | 0.000e+00 | 1.021e-19 | 0.000e+00 | 1.661e-12 | True | 1,1,1 |
| 3 | 0.20 | 20261005 | finufft | 11.455 | 5.735e-20 | 0.000e+00 | 8.258e-20 | 0.000e+00 | 1.976e-12 | True | 1,1,1 |
| 3 | 0.20 | 20261006 | finufft | 11.692 | 4.709e-20 | 0.000e+00 | 8.401e-20 | 0.000e+00 | 1.796e-12 | True | 1,1,1 |
| 3 | 0.20 | 20261007 | finufft | 13.439 | 3.369e-20 | 0.000e+00 | 1.814e-19 | 0.000e+00 | 1.566e-12 | True | 1,1,1 |
| 3 | 0.20 | 20261008 | finufft | 11.772 | 4.280e-20 | 0.000e+00 | 7.085e-20 | 0.000e+00 | 1.789e-12 | True | 1,1,1 |
| 3 | 0.20 | 20261009 | finufft | 12.139 | 6.837e-20 | 0.000e+00 | 3.007e-19 | 0.000e+00 | 1.958e-12 | True | 1,1,1 |
| 3 | 0.20 | 20261010 | finufft | 12.192 | 2.975e-20 | 0.000e+00 | 1.310e-19 | 0.000e+00 | 1.920e-12 | True | 1,1,1 |
| 3 | 0.30 | 20261001 | finufft | 18.216 | 3.561e-20 | 0.000e+00 | 8.793e-20 | 0.000e+00 | 2.023e-12 | True | 1,1,1 |
| 3 | 0.30 | 20261002 | finufft | 20.625 | 5.876e-20 | 0.000e+00 | 2.257e-19 | 0.000e+00 | 1.723e-12 | True | 1,1,1 |
| 3 | 0.30 | 20261003 | finufft | 19.155 | 4.397e-20 | 0.000e+00 | 1.605e-19 | 0.000e+00 | 1.800e-12 | True | 1,1,1 |
| 3 | 0.30 | 20261004 | finufft | 20.737 | 6.508e-20 | 0.000e+00 | 2.618e-19 | 0.000e+00 | 1.658e-12 | True | 1,1,1 |
| 3 | 0.30 | 20261005 | finufft | 18.136 | 7.319e-20 | 0.000e+00 | 1.752e-19 | 0.000e+00 | 2.024e-12 | True | 1,1,1 |
| 3 | 0.30 | 20261006 | finufft | 19.605 | 3.598e-20 | 0.000e+00 | 8.137e-20 | 0.000e+00 | 1.452e-12 | True | 1,1,1 |
| 3 | 0.30 | 20261007 | finufft | 18.905 | 3.344e-20 | 0.000e+00 | 8.056e-20 | 0.000e+00 | 1.619e-12 | True | 1,1,1 |
| 3 | 0.30 | 20261008 | finufft | 18.759 | 3.759e-20 | 0.000e+00 | 6.038e-20 | 0.000e+00 | 1.725e-12 | True | 1,1,1 |
| 3 | 0.30 | 20261009 | finufft | 19.777 | 6.405e-20 | 0.000e+00 | 2.049e-19 | 0.000e+00 | 2.074e-12 | True | 1,1,1 |
| 3 | 0.30 | 20261010 | finufft | 19.948 | 5.122e-20 | 0.000e+00 | 1.692e-19 | 0.000e+00 | 1.358e-12 | True | 1,1,1 |
| 3 | 0.40 | 20261001 | finufft | 46.505 | 2.754e-20 | 0.000e+00 | 9.382e-20 | 0.000e+00 | 1.682e-12 | True | 5,1,1 |
| 3 | 0.40 | 20261002 | finufft | 44.862 | 3.916e-20 | 0.000e+00 | 1.173e-19 | 0.000e+00 | 1.436e-12 | True | 5,1,1 |
| 3 | 0.40 | 20261003 | finufft | 45.927 | 6.660e-20 | 0.000e+00 | 1.465e-19 | 0.000e+00 | 1.569e-12 | True | 5,1,1 |

## Interpretation and scope

Failure to find a candidate within this finite random-start budget does not prove infeasibility. Validation also tests five translated 20x20 periodic lattices at all seven shell selections; those constructed candidates are separate from random-start success rates. No chi=0.60, adaptive refinement, added seeds, or continuation stage is included.

Raw receipts, per-run initial/final coordinates and species labels, complete mode lists, direct-reference structure factors, optimizer codes, and backend-use counters are retained. state.json records the current job and any blocking condition. The service is a finite Python process, not an agent coordinator.
