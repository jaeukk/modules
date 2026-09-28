# One-initial-state diagnosis: energy contributions, component drift, and stopping tolerance

## Main result

The requested controls do not resolve the vanilla–FINUFFT endpoint discrepancy for seed 20260928. Increasing the allowed residual energy has almost no effect on the endpoint separation, and independent component drift accounts for very little of it. The traces locate the separation during the early coupled collective-coordinate/contact relaxation, far above ground-state energies.

This task does not change the optimizer, physical objective, production code, or FINUFFT accuracy. Only one initial state was used: N=2000, two dimensions, five equal components, chi per component=80/798, total packing fraction=0.2, four OpenMP threads, one internal FINUFFT thread, FINUFFT epsilon=1e-12. The existing fixed-order contact sum is retained.

## 1. Full trajectories and separate energies

The existing matched baseline trajectories were reused and fully replayed: 1141 vanilla and 1150 FINUFFT objective evaluations. Every replay reproduced the saved native total energy exactly. At every configuration we recorded native collective-coordinate energy and independently computed contact energy using all particle pairs; their sum was checked against the native total. Independent direct collective-coordinate energy, overlap, and force checks were also recorded.

Every particle position is saved for every callback in the original binary traces and portable compressed NPZ files. Callback indices include line-search trial points: **they are not accepted iterations or physical time**. Comparing the same callback number after the paths diverge does not establish identical optimizer stages.

| Callback | Vanilla CC energy | FINUFFT CC energy | Vanilla contact energy | FINUFFT contact energy | Mean position difference, Cartesian |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 15.6911552823215 | 15.6911552823219 | 121.534198072801 | 121.534198072801 | 0 |
| 40 | 4.957054480202 | 4.957054480197 | 0.141686362416 | 0.141686362502 | 5.69e-12 |
| 50 | 3.935788708708 | 3.935788705128 | 0.159221858694 | 0.159221880017 | 5.49e-10 |
| 57 | 3.290777633641 | 3.290776786675 | 0.040595962887 | 0.040597277233 | 1.63e-7 |
| 68 | 2.741368515747 | 2.741291712812 | 0.095769032165 | 0.095997142395 | 1.83e-5 |
| 100 | 1.656976787500 | 1.734465809989 | 0.141069797620 | 0.017243453928 | 0.026834 |
| 800 | 0.001433669070 | 0.002031826382 | 4.24e-6 | 3.33e-6 | 0.420420 |

At the identical initial coordinates, contact energies agree; the approximately 4.26e-13 total-energy difference is in the collective-coordinate contribution. After positions begin to separate, contact-energy differences also grow. For example, at callback 50 the contact-energy difference is approximately 2.13e-8, versus a collective-coordinate difference of -3.58e-9.

Contact identities are identical through callback 73. The first differing contact pair appears at callback 74, when the mean Cartesian displacement is already 7.82e-5. By callback 100, 77 pairs differ. Thus the first contact-set difference does not precede the initial positional separation in these callback records. This does **not** rule out contact-force curvature or contact distances within the unchanged set as contributors to amplification.

The energy traces show intermittent contact-energy excursions while collective-coordinate energy is reduced. These include optimizer trial evaluations and should not be interpreted as increases along an accepted-step sequence. The results support amplification during coupled relaxation, not an explanation based only on final near-zero energy arithmetic. They do not yet identify the exact accepted-step/line-search decision that separates the paths or prove that contact repulsion alone causes the sensitivity.

![Energy, component displacement, and contact traces](results/trajectories.png)

[Early-energy detail, first 200 callbacks](results/energies_first_200.png)

## 2. Increasing the ground-state stopping tolerance

Here “increase tolerance” means a larger positive total-energy stopping target. We used the existing `RelaxStructure_NLOPT_Emin` routine, keeping L-BFGS, history 100, the three coordinate scales, and FINUFFT epsilon unchanged. These are eight actual early-stop optimizations, not post-hoc changes to the ground-state acceptance label.

For every run, the recorded callback trajectory exactly matches the corresponding prefix of the original baseline trace. All starting coordinates are byte-identical. Actual backend counters confirm FINUFFT use without optimization fallback.

| Total-energy stop target | Vanilla evaluations | FINUFFT evaluations | Vanilla seconds | FINUFFT seconds | Mean endpoint separation, Cartesian |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1e-4 | 929 | 948 | 5.427 | 5.762 | 0.431369213 |
| 1e-8 | 1072 | 1092 | 6.654 | 5.813 | 0.430885934 |
| 1e-11 | 1101 | 1107 | 6.374 | 6.147 | 0.430885985 |
| 1e-14 | 1114 | 1121 | 6.413 | 6.256 | 0.430885974 |
| Original baseline, final energy about 1e-19 | 1141 | 1150 | — | — | 0.430885974 |

Timings include position/gradient trace writing and are not isolated performance benchmarks. The 1e-4 and 1e-8 states fail the original stricter ground-state acceptance criteria; they are intentionally looser stopping controls, not newly certified ground states. The 1e-11 and 1e-14 states satisfy those criteria.

Interpretation: the trajectories have already separated substantially before any of these energy targets is reached. Stopping earlier in this range cannot undo that separation.

## 3. Component-resolved displacements

Particles retain their identities and species labels. Differences are periodic and measured in Cartesian units. Each component's translation is the mean wrapped displacement of its particles; removing it isolates internal rearrangement. The small observed displacements relative to the box support this local translation estimate without ambiguous winding.

| Component | Raw mean displacement | Mean after its own drift removal | Component translation vector, Cartesian |
| --- | ---: | ---: | --- |
| 0 | 0.418135 | 0.417854 | (0.017029, -0.004535) |
| 1 | 0.421785 | 0.421561 | (-0.020420, 0.011590) |
| 2 | 0.386962 | 0.386446 | (-0.025506, 0.008976) |
| 3 | 0.452205 | 0.452033 | (0.029617, 0.008973) |
| 4 | 0.475343 | 0.475103 | (-0.000720, -0.025004) |

The equal-component average decreases from 0.430886 to 0.430599 after separate drift removal, only about 0.067%. Relative component drift exists but does not explain the large endpoint mismatch.

## Artifacts and reproduction

- [Machine-readable summary](results/summary.json), [component endpoints](results/endpoint_components.json)
- [Vanilla energy history](results/direct_energies.csv), [FINUFFT energy history](results/finufft_energies.csv)
- [Position differences and contact identities by callback](results/first_pass_displacements_contacts.csv)
- [Vanilla full particle histories](results/direct_positions.npz), [FINUFFT full particle histories](results/finufft_positions.npz)
- Each stopping run retains initial/final positions, a full callback trace, direct structure factors, and its command/result receipt.
- `diagnose.cpp` provides `decompose` and `optimize` modes. `build.py` compiles against hash-checked prior core objects; `run.py` executes and analyzes this one-seed study. `run.py` refuses to overwrite the existing `results` directory.
- [Evidence manifest](results/manifest.json) records settings, source/binary hashes, and hashes of the reused baseline inputs. Current build commands are saved by `build.py`.

Production files and the previous multiagent audit artifacts are unchanged. No alternative optimizer or additional random initial state was tested in this task. The cross-backend endpoint agreement issue remains unresolved.
