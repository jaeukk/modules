# Statistics — n-body correlation functions $g_n$ of 1D point patterns

Computes the ensemble-averaged **n-body correlation function** $g_n$ ($n = 2, 3, 4$)
and the **$\tau_n$ order metric** of one-dimensional point patterns stored in a
`ConfigurationPack`, plus hole and overlap statistics. Migrated from the legacy
`new_c++/Statistics/` consumer; all configurations of an ensemble are assumed to have
the same box length. Validation tools are kept apart in [tests/](tests/).

## Layout

| file | contents |
|------|----------|
| [gn.h](gn.h) / [gn.cpp](gn.cpp) | `gn`: binned $g_n$ estimator with checkpointing (resumable), $\tau_n$ and long-range $\tau_n$ from the stored $g_n$. |
| [tau_n.h](tau_n.h) / [tau_n.cpp](tau_n.cpp) | `tau_n`: $\tau_n$ without storing $g_n$ (global-grid estimator; legacy algorithm, not optimized). |
| [main.cpp](main.cpp) | `gn.out` CLI: modes `gn_1D`, `taun_1D`, `hole`, `overlap`, `test`. |
| [tests/](tests/) | Not part of the app: test-input generator `gen.out` and the exact-$g_n$ check `check_url_gn.py` (see [tests/README.md](tests/README.md)). |

## Build

```sh
make            # -> gn.out (builds core objects first)
make clean
make -C tests   # validation tools (adds models)
```

Depends on `core`; GSL, OpenMP.

## $g_n$ estimator (`gn_1D`)

For every particle $x_1$ the relative positions $x_{j1} = x_j - x_1$ of all neighbors
(periodic images included) are histogrammed into bins of width $\Delta x$ centered at
$i\,\Delta x$, $|i| \le n_h = \lfloor R_\max/\Delta x \rfloor$. With $N_i$ the occupancy
of bin $i$ ($x_1$ itself sits in bin 0),

$$
g_n(x_{21},\dots,x_{n1}) \approx \frac{1}{N\rho^{n-1}\Delta x^{n-1}}
\sum_{x_1}\prod_{k=2}^{n}\Big(N_{i_k} - \#\{j<k : i_j = i_k\}\Big),
$$

i.e. the bin average of $g_n$ over a cube of side $\Delta x$ (ordered tuples of
*distinct* particles). Only $x_{21} \ge 0$ is stored; $g_2$ is symmetrized over
$\pm x_{21}$. Memory $\propto (R_\max/\Delta x)^{n-1}$; time
$\propto N_c\,N\,(2\rho R_\max)^{n-1}$. The ensemble mean and its standard error
(over configurations) are written.

### Input (stdin), `./gn.out gn_1D [hours] [n]`

```
<n>                 # 2, 3, or 4
<checkpoint prefix> # progress is saved to <prefix>.g<n>; an existing file is resumed
<Rmax> <dR>         # ignored when resuming
<num threads>
<ConfigPack prefix>
<num configs>
<output prefix>
```

`hours` is the wall-clock limit (default 10); the checkpoint is written when it is near,
and rerunning the same input continues. A third argument `n` skips writing the $g_n$
text files. Outputs: `<out>_g2.txt` / `<out>_g3.txt` (columns $x_{21}$, [$x_{31}$,] $g_n$,
error) or one `<out>_g4_x21_<x21>.txt` per $x_{21}$ bin (columns $x_{31}$, $x_{41}$, $g_4$,
error); `<out>_tau<n>` and `<out>_tau<n>_exc` ($\tau_n(L)$ and the variant excluding
pair separations below `Rmin` = 5).

### Changes from the legacy code

| change | legacy behavior |
|--------|-----------------|
| Neighbors searched to $(n_h + \tfrac12)\Delta x$; every bin fully sampled | outermost bins $|i| = n_h$ half-filled (biased ≈ ×½); last $g_2$ bin always 0 |
| Normalized by $\rho^{n-1}$ | assumed $\rho = 1$ |
| Bins centered at $i\,\Delta x$ for any $R_\max/\Delta x$ | non-integer $R_\max/\Delta x$ aborted (`std::out_of_range`) |
| `gn_1D <hours> n` honored | tested `argv[2]`, so `n` was ignored |
| Long-range $\tau_n$: no `result[-1]` read | out-of-bounds read when `Rmin` ≤ 2Δx |
| Invalid `n` → error exit | half-constructed object, then crash |

Interior bins are **unchanged** (identical printed values vs. the legacy build for
$n=2,3,4$). The checkpoint layout is unchanged, so legacy `.g<n>`
files resume (their edge bins keep the legacy bias).

### Performance

The legacy kernel opened one OpenMP region per particle ($n=4$) and swept the whole
$(R_\max/\Delta x)^{n-1}$ array every configuration. Now: sorted-position window scan
instead of the generic cell list; one parallel region per configuration (none below
~$10^5$ operations), dynamically scheduled over $x_{21}$ slabs that each thread owns
(no atomics); only touched bins are accumulated. For $n=4$ the counts are symmetric
under $x_{31} \leftrightarrow x_{41}$, so only $x_{31} \le x_{41}$ is accumulated (halves the
work; outputs byte-identical). Results do not depend on the thread count (byte-identical
checkpoints for 1 and 4 threads).

Cloaked URL, $N=1000$, 4 threads, on a loaded 8-core WSL machine (absolute times are
noisy):

| case | legacy | now |
|------|--------|-----|
| $n=2$, $R_\max=3$, $\Delta x=0.05$, 1000 configs | 12 s | 0.16 s |
| $n=3$, $R_\max=3$, $\Delta x=0.1$, 1000 configs | 12 s | 0.29 s |
| $n=4$, $R_\max=2$, $\Delta x=0.1$, 1000 configs | 697 s (8 threads) | 0.28 s |
| $n=4$, $R_\max=10$, $\Delta x=0.1$ (4.1 M bins), per config | 50 ms | 8.6 ms → 4.2 ms with the $x_{31} \leftrightarrow x_{41}$ symmetry (idle machine) |

Cost per configuration scales with the number of occupied-bin tuples,
$\approx N m^{n-1}/2$ ($m$ = occupied bins around a particle, $\le$ neighbors + 1), not with
the number of bins: per configuration at most that many bins are touched (48 % of the
4.1 M bins at $\Delta x = 0.1$, 9 % at $\Delta x = 0.05$), and the $\sum v,\ \sum v^2$ update
costs about half of the kernel.

## Validation

Checked against the exact $g_n$ of the 1D cloaked uniformly randomized lattice
($a = 1, 2$; $n = 2, 3, 4$; $\rho = 1$ and $1/2$): all PASS. Method, commands, and
results: [tests/README.md](tests/README.md).

## Other modes

- `taun_1D [hours]` — stdin: `n ConfigPack Nc threads out_prefix ckpt_prefix Lmax dx`.
- `hole` — hole probability $H_V$ and $E_V$; stdin: `ConfigPack begin end threads dR samples prefix`.
- `overlap` — overlap function $Q(t;a)$; stdin: `ConfigPack begin end threads sigma_max dsigma prefix`.
- `test` — index ↔ multi-index round trip.
