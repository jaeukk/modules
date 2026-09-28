# Five-component FINUFFT internal-thread control

N=2000, five equal components, chi_i=80/798, phi_total=0.2. Fixed-order contact sum. OMP_NUM_THREADS=OMP_THREAD_LIMIT=4; CPU affinity 0,2,4,6. Only the options passed to finufft2d1 are changed: opts.nthreads=1. The link wrapper records every override and status. Original object files, optimizer, tolerance 1e-12, and other options are unchanged.

Three seeds, two separate process repeats per backend. Initial files and mode files match the preceding four-thread FINUFFT experiment byte for byte. Times exclude independent endpoint validation. Displacements are Cartesian, fixed-label, periodic, and translation-aligned.

| Seed | Vanilla seconds r0/r1 | FINUFFT seconds r0/r1 | Vanilla repeat mean | FINUFFT repeat mean | Cross-backend mean r0/r1 |
|---|---|---|---:|---:|---|
| 20260928 | 7.087/7.015 | 6.260/6.325 | 0.000000e+00 | 0.000000e+00 | 4.308860e-01/4.308860e-01 |
| 20260929 | 7.275/7.289 | 7.138/6.487 | 0.000000e+00 | 0.000000e+00 | 4.106835e-01/4.106835e-01 |
| 20260930 | 6.549/6.426 | 5.532/5.547 | 0.000000e+00 | 0.000000e+00 | 2.863720e-01/2.863720e-01 |

All ground-state checks passed: True. All optimizer codes successful: True.
Repeat endpoints byte-identical: {'direct': True, 'finufft': True}.

| Backend | Median optimization seconds | Minimum final energy | Maximum final energy |
|---|---:|---:|---:|
| direct | 7.050934 | 2.272589e-20 | 9.327793e-20 |
| finufft | 6.292521 | 2.448165e-20 | 1.383751e-19 |

Same-backend reproducibility and agreement between backends are distinct tests. This finite control does not trace the amplification mechanism or guarantee repeatability on other machines or configurations. No production thread default was changed; the broader sweep remains stopped.
