# Statistics/tests — validation of `gn.out`

Tools to check the $g_n$ estimator of `../gn.out` against an exactly solvable model.
They are not part of the app and are built separately.

| file | contents |
|------|----------|
| [gen.cpp](gen.cpp) | `gen.out`: test inputs — perturbed lattices (`PL`) and Poisson patterns (`Poisson`), written as a ConfigPack. Migrated from the legacy `new_c++/Statistics/gen.cpp`. |
| [check_url_gn.py](check_url_gn.py) | Compares `gn_1D` output with the **exact** $g_n$ of the uniformly randomized lattice (URL), or two runs with each other. |

## Build

```sh
make            # -> gen.out (builds core / models objects first)
make clean
```

Depends on `core` + `models`; GSL, OpenMP. `check_url_gn.py` needs only NumPy (`uv run`
resolves it from the inline metadata).

## Test inputs (`gen.out`)

- `gen.out PL` — stdin: `dim uniform <R> L Nc prefix` (seed 0). Each site of the
  integer/square/cubic lattice with `L` cells per side moves by a vector uniform in a ball
  of radius `R`; in 1D this is the URL with $a = 2R$, which is cloaked for integer $a$.
- `gen.out Poisson` — stdin: `rho mode(0 canonical/1 grand canonical) dim volume Nc prefix`.

The reference ensembles (1D, $N = L = 1000$, 1000 configurations):

```sh
printf "1\nuniform\n0.5\n1000\n1000\nURL_a1_N1000\n" | ./gen.out PL   # a = 1
printf "1\nuniform\n1.0\n1000\n1000\nURL_a2_N1000\n" | ./gen.out PL   # a = 2
```

## Exact check (`check_url_gn.py`)

For the 1D URL ($u_k \sim U[-a/2, a/2)$ i.i.d., integer $a$ ⇒ cloaked; Klatt, Kim &
Torquato, PRE **101**, 032118 (2020)),

$$
g_n(0, r_2, \dots, r_n) = a^{-(n-1)}\,\mathbb{E}_u\Big[\#\{\text{injective }(k_2..k_n):
k_j \ne 0,\ |r_j + u - k_j| < a/2\}\Big],\qquad u \sim U[-a/2, a/2),
$$

integrated exactly between its breakpoints in $u$ (for $a=1$: $g_2 = \min(|r|,1)$).
The script bin-averages it and reports $z = (g_n - g_n^\text{exact})/\delta g_n$:

```sh
printf "3\nck\n3\n0.1\n4\nURL_a1_N1000\n1000\nout\n" | ../gn.out gn_1D 1
uv run check_url_gn.py out --n 3 --a 1 --sub 6     # vs exact: PASS/FAIL (--sub: reference sub-samples per bin)

# Compare two runs, e.g. the same input with another gn.out build (output prefix "ref"):
printf "3\nck_ref\n3\n0.1\n4\nURL_a1_N1000\n1000\nref\n" | ../gn.out gn_1D 1
uv run check_url_gn.py out --n 3 --compare ref     # add --interior to skip the outermost bins
```

A checkpoint prefix (`ck`, `ck_ref`) that already exists is resumed, not recomputed; use a
fresh prefix for a new run.

The exact test cannot see errors confined to bins where $g_n < 0.05$ (bins without samples
are tolerated there), e.g. a dropped $x_{31} = x_{41}$ diagonal for $a = 1$. After changing the
kernel, also `--compare` against a run of the previous build.

## Status

1000 configurations, $N = 1000$, $a = 1, 2$, $n = 2, 3, 4$: all PASS, rms $z$ = 0.96–1.04
over all bins. Pooled bias $\langle g_n - g_n^\text{exact}\rangle$ is within 1.5 standard
errors in every case; $\tau_n(L_\max)$ agrees with the exact value to 0.1 %. The
$\rho^{n-1}$ normalization was checked on a spacing-2 URL ($\rho = 1/2$): far-field
$g_2 = 1.0001(12)$, $g_3 = 0.9992(8)$ (the legacy $\rho = 1$ normalization would give 0.50
and 0.25). The exact formula was cross-checked against brute-force Monte Carlo, and the
$n=4$ estimator against direct enumeration of ordered quadruples (agreement to
$2\times10^{-16}$). Bins where $g_n$ is tiny have few, overdispersed counts (one cluster of
$n$ particles fills several bins), so empty bins there are expected.
