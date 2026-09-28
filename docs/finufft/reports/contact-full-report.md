# Contact-repulsion and fixed-order summation study

Three matched seeds, two independent process repeats per backend. All runs use phi_total=0.2, equal radii sqrt(0.2/pi), density=1, u0=e0=1, and the unchanged three-pass L-BFGS routine (2000 evaluations/pass, 600 seconds/configuration). One-component N500 uses the previous 299 modes (chi=299/998); five-component N2000 uses 80 complete-shell modes per species (chi=80/798).

Times are optimization wall times. Each endpoint is independently checked by direct Fourier sums and all-pair contact forces. Displacements use fixed labels, periodic wrapping, and optimal translation. Acceptance is E<=1e-11, maxS<=1e-10, max force component<=1e-7, max relative overlap<=1e-6.

## baseline

| Case | Seed | Vanilla seconds (r0/r1) | FINUFFT seconds (r0/r1) | Vanilla repeat mean | FINUFFT repeat mean | Cross-backend mean (r0/r1) | All accepted |
|---|---|---|---|---:|---:|---|---|
| one_2t | 20260928 | 1.324/1.296 | 0.923/0.927 | 0.000000e+00 | 5.305981e-13 | 1.566474e-11/1.526783e-11 | True |
| one_2t | 20260929 | 1.264/1.264 | 0.908/0.938 | 0.000000e+00 | 2.384501e-13 | 5.768844e-11/5.772638e-11 | True |
| one_2t | 20260930 | 1.292/1.256 | 0.918/0.914 | 0.000000e+00 | 5.613782e-13 | 2.624895e-11/2.623255e-11 | True |
| one_4t | 20260928 | 0.713/0.749 | 0.681/0.688 | 0.000000e+00 | 1.276130e-12 | 1.555872e-11/1.581433e-11 | True |
| one_4t | 20260929 | 0.721/0.753 | 0.702/0.737 | 0.000000e+00 | 2.432929e-12 | 5.630546e-11/5.428698e-11 | True |
| one_4t | 20260930 | 0.807/0.698 | 0.740/0.707 | 0.000000e+00 | 2.794406e-13 | 2.635596e-11/2.613379e-11 | True |
| five_4t | 20260928 | 6.538/6.589 | 10.091/9.865 | 1.829184e-02 | 3.824413e-01 | 4.235938e-01/4.077442e-01 | True |
| five_4t | 20260929 | 7.304/7.172 | 10.323/10.672 | 2.603880e-02 | 3.707791e-01 | 3.927150e-01/3.974823e-01 | True |
| five_4t | 20260930 | 6.476/6.363 | 9.896/10.765 | 4.357369e-02 | 2.813897e-01 | 2.893144e-01/2.797577e-01 | True |
## fixed

| Case | Seed | Vanilla seconds (r0/r1) | FINUFFT seconds (r0/r1) | Vanilla repeat mean | FINUFFT repeat mean | Cross-backend mean (r0/r1) | All accepted |
|---|---|---|---|---:|---:|---|---|
| one_2t | 20260928 | 1.360/1.335 | 0.936/0.873 | 0.000000e+00 | 1.778119e-13 | 1.545069e-11/1.545872e-11 | True |
| one_2t | 20260929 | 1.291/1.340 | 0.911/0.976 | 0.000000e+00 | 9.733354e-13 | 5.766966e-11/5.776289e-11 | True |
| one_2t | 20260930 | 1.270/1.215 | 0.910/0.923 | 0.000000e+00 | 1.881961e-12 | 2.696342e-11/2.675999e-11 | True |
| one_4t | 20260928 | 0.725/0.710 | 0.690/0.672 | 0.000000e+00 | 1.424247e-12 | 1.600835e-11/1.543678e-11 | True |
| one_4t | 20260929 | 0.726/0.735 | 0.683/0.684 | 0.000000e+00 | 2.034196e-13 | 5.426136e-11/5.427905e-11 | True |
| one_4t | 20260930 | 0.686/0.684 | 0.650/0.662 | 0.000000e+00 | 2.562093e-13 | 2.642786e-11/2.649684e-11 | True |
| five_4t | 20260928 | 7.490/6.978 | 10.308/9.912 | 0.000000e+00 | 4.238511e-01 | 4.289947e-01/4.208280e-01 | True |
| five_4t | 20260929 | 7.023/8.137 | 10.469/10.160 | 0.000000e+00 | 4.014114e-01 | 3.794657e-01/3.987041e-01 | True |
| five_4t | 20260930 | 6.678/7.262 | 10.324/9.605 | 0.000000e+00 | 2.924702e-01 | 2.954983e-01/2.957894e-01 | True |

Per-run JSON files retain energies, forces, maximum overlap, optimizer statuses and evaluation counts, transform counters, timing, and input hashes. Probe JSON files retain 100 repeated fixed-coordinate evaluations plus independent-force and finite-difference checks. The baseline and fixed executables and their build commands/source hashes are preserved separately.
