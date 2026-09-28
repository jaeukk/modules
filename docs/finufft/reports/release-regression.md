# FINUFFT nthreads=1 repeatability regression

Result: **PASS** (101 required checks, 0 failed).

N=2000, five equal components, phi_total=0.2, chi_i=80/798, OMP_NUM_THREADS=OMP_THREAD_LIMIT=4, passive waiting, CPU affinity [0, 2, 4, 6]. Production defaults only: the FINUFFT option is never set by the harness, so these runs exercise the default nthreads=1.

Finite regression evidence on this machine, build, and settings only; not a general guarantee. Cross-backend displacements are descriptive, not acceptance criteria.

| Run | Actual backend | Optimization s | Process s | Evaluations | Status codes | Accepted | FINUFFT calls (harness/observed) |
|---|---|---:|---:|---:|---|---|---|
| 20260928_direct_0 | direct | 6.729 | 6.808 | 1141 | 1,1,1 | True | 0/0 |
| 20260928_finufft_0 | finufft | 6.664 | 6.750 | 1150 | 1,1,1 | True | 5765/5770 |
| 20260928_finufft_1 | finufft | 6.317 | 6.403 | 1150 | 1,1,1 | True | 5765/5770 |
| 20260928_direct_1 | direct | 6.970 | 7.046 | 1141 | 1,1,1 | True | 0/0 |
| 20260929_finufft_0 | finufft | 6.206 | 6.289 | 1146 | 1,1,1 | True | 5745/5750 |
| 20260929_direct_0 | direct | 6.894 | 6.959 | 1135 | 1,1,1 | True | 0/0 |
| 20260929_direct_1 | direct | 6.859 | 6.903 | 1135 | 1,1,1 | True | 0/0 |
| 20260929_finufft_1 | finufft | 5.997 | 6.053 | 1146 | 1,1,1 | True | 5745/5750 |
| 20260930_direct_0 | direct | 6.059 | 6.101 | 1046 | 1,1,1 | True | 0/0 |
| 20260930_finufft_0 | finufft | 5.472 | 5.546 | 988 | 1,1,1 | True | 4955/4960 |
| 20260930_finufft_1 | finufft | 5.433 | 5.503 | 988 | 1,1,1 | True | 4955/4960 |
| 20260930_direct_1 | direct | 6.312 | 6.398 | 1046 | 1,1,1 | True | 0/0 |

| Seed | Direct repeats identical | FINUFFT repeats identical | Cross-backend mean displacement r0/r1 |
|---|---|---|---|
| 20260928 | True | True | 4.309e-01/4.309e-01 |
| 20260929 | True | True | 4.107e-01/4.107e-01 |
| 20260930 | True | True | 2.864e-01/2.864e-01 |

direct: median optimization 6.794 s, range 6.059-6.970 s.

finufft: median optimization 6.102 s, range 5.433-6.664 s.
