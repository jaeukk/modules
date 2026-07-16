# OrderMetric_PeriodicNetwork

Two independent programs for two-phase-media descriptors of nonspherical
(polygonal / polyhedral) systems:

1. **`NetworkBV.out`** — the **local volume-fraction variance** $\sigma_V^2(R)$ and
   the **hyperuniformity order metric** $\overline{B}_V$ of a periodic network
   (a two-phase medium built from a predefined 2D/3D lattice of polygonal or
   spherical struts/particles).
2. **`GetChivDispersions.out`** — the **spectral density** $\tilde{\chi}_V(k)$ of a packing
   of **polygonal (2D)** and **polyhedral (3D)** particles.

Both are migrated from the legacy `NonsphericalParticles/` project and depend on
the `nonspherical` module (which extends the sphere-packing descriptors to
arbitrary shapes; see `modules/nonspherical/README.md`).

## Provenance

Program **(1)**, applied to the predefined periodic networks, was used to estimate
the data published in

> J. Kim & S. Torquato, *Characterizing the hyperuniformity of ordered and
> disordered two-phase media*, **Phys. Rev. E 103**, 012123 (2021).

There $\overline{B}_V$ (the coefficient of the small-$R$–independent, large-$R$
$R^{-(d+1)}$ decay of $\sigma_V^2(R)$) is the order metric that ranks the
hyperuniformity of the two-phase media.

## Build

```sh
make            # -> GetChivDispersions.out AND NetworkBV.out  (builds core + nonspherical + pairstat + models objects first)
make clean
```

`-DUSE_SHAPES` is set automatically by the makefile.

---

## Program 1 — `NetworkBV.out` (network $\sigma_V^2(R)$ and order metric $\overline{B}_V$)

Interactive: reads answers on **stdin**, echoes each prompt to stdout. It builds a
periodic network at a sequence of volume fractions and computes $\overline{B}_V$
either in **Fourier space** (a reciprocal-lattice sum; required for polygonal
particles) or in **direct space** (Monte-Carlo local volume-fraction variance;
spheres only).

**Predefined networks** (the `Network name` prompt):

- 2D: `SquareNetwork`, `OctagonalNetwork`, `RhombicNetwork`, `TriangularNetwork`,
  `HoneycombNetwork`, `KagomeNetwork`
- 3D: `CubicNetwork`, `FCCNetwork`, `BCCNetwork`, `HCPNetwork` (and Wigner–Seitz
  variants `FCCNetwork_WS` / `BCCNetwork_WS` / `HCPNetwork_WS`)

**Prompt sequence** (in order):

| prompt | value |
|--------|-------|
| `Particle Type (sphere/polygon)` | `polygon` (2D/3D shape network) or `sphere` |
| `Do compute Bv in the direct space? (sphere only)` | `0` (Fourier) or `1` (direct; sphere only) |
| *Fourier:* `Maximal / Minimal wavenumber`, `Change in wavenumber` | `kmax kmin dk` |
| *(or Direct, sphere only:* `Maximal radius`, `dR`, `Minimal radius`, `Interval of R`*)* | `Rmax dR Rmin intv_R` |
| `The number of threads` | e.g. `4` |
| `Initial / Final / Change in volume fraction` | `phi_init phi_fin d_phi` |
| `Network name` | e.g. `SquareNetwork` |
| `Lengthscale (D/ rho/ a)` then its value | `a` + lattice constant (or `D`+nn-distance, `rho`+number density) |

It then reports $\overline{B}_V$ over the requested $\phi$ range and asks whether
to process another cell (answer `n` to stop).

**Example** (2D square network of polygonal struts, Fourier space):

```sh
./NetworkBV.out <<'EOF'
polygon
0
1000
100
50
4
0.05
0.50
0.05
SquareNetwork
a
1.0
n
EOF
```

**Outputs:** `<network>_Bv_Fourier.txt` (columns: $\phi$, $\overline{B}_V$) — or
`<network>_Bv_Direct.txt` plus `<network>_runningav_<phi>` files for the
direct-space (sphere) mode.

---

## Program 2 — `GetChivDispersions.out` (spectral density $\chi_V(k)$ of polygon/polyhedron packings)

Reads an ensemble of nonspherical **dispersions** from disk (`.Dispersion` files,
via `ReadDispersion`) and computes the orientationally-averaged spectral density
$\chi_V(k)$, using each particle's shape form factor.

**Usage:**

```sh
./GetChivDispersions.out <load_prefix> <num_configs> <save_prefix> <circular_kmax> <linear_kmax> <dk>
```

| arg | meaning |
|-----|---------|
| `load_prefix` | printf-style prefix expanded with the config index to locate each `.Dispersion` file (e.g. `config_%04d`) |
| `num_configs` | number of configurations to average over |
| `save_prefix` | output prefix for the $\chi_V(k)$ table |
| `circular_kmax` | all reciprocal-lattice vectors with $\lvert k\rvert$ below this are summed |
| `linear_kmax` | multiples of those k-vectors are taken up to this $\lvert k\rvert$ |
| `dk` | binning width in $\lvert k\rvert$ |

**Example** (100 dispersions named `poly_0000`, `poly_0001`, …):

```sh
./GetChivDispersions.out poly_%04d 100 chiv_poly 20 40 0.1
# -> chiv_poly : columns (k, chi_V(k), error, ...)
```

Two convenience modes:
- `./GetChivDispersions.out <prefix>` — read a **single** dispersion.
- `./GetChivDispersions.out` (no args) — run the interactive **form-factor test**
  (analytic vs. numerical form factors of a triangle / square / cube).

**Output:** a $\chi_V(k)$ table written to `<save_prefix>`.

---

## Dependencies & notes

- **Modules:** `core`, `nonspherical` (both programs); `NetworkBV.out` additionally
  links `pairstat` (spherical `SpectralDensity`, `MC_vfv`) and `models`
  (`LatticeGen`). All include/library paths come from `../../config.mk`.
- **External:** Boost (header-only), OpenMP. No CGAL/voro++/FFTW — runnable under
  default flags.
- The spectral-density kernel `SpectralDensity_Dispersion` lives in the
  **`nonspherical` module** (linked here), not in this folder.
- The folder name is historical: `main.cpp` is the $\chi_V(k)$ driver above, **not**
  an order-metric driver. The legacy order-metric `main()` was `Test_ParticleShapes.cpp`,
  which was not migrated.
