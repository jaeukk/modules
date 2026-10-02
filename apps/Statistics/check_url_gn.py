# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy"]
# ///
"""Check gn.out output for a 1D uniformly randomized lattice (URL) against the exact g_n.

URL (Klatt, Kim & Torquato, PRE 101, 032118 (2020)): every site k of the integer
lattice is displaced by u_k ~ U[-a/2, a/2), i.i.d.  For integer a the URL is
cloaked, and conditioning on a particle at the origin with displacement u gives

    g_n(0, r_2, ..., r_n) = a^{-(n-1)} E_u[ #{injective (k_2..k_n) : k_j != 0,
                                              |r_j + u - k_j| < a/2} ],

which is piecewise constant in u and is integrated exactly between breakpoints.
For a = 1 this reduces to g_2(r) = min(|r|, 1).

Usage:
    uv run check_url_gn.py PREFIX --n 3 --a 1            # estimate vs exact
    uv run check_url_gn.py PREFIX --n 3 --compare OTHER  # two runs vs each other
"""

from __future__ import annotations

import argparse
import itertools
import re
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

ZERO_ERR_TOL = 0.05  # largest |g_n - ref| tolerated in a bin that has no samples


def gn_url_exact(r: NDArray, a: int) -> NDArray:
    """Exact g_n of the 1D URL with integer perturbation strength a.

    Parameters
    ----------
    r : (M, n-1) ndarray
        Relative positions x_j - x_1 (j = 2..n) [lattice constant].
    a : int
        Side length of the displacement interval [-a/2, a/2) [lattice constant].

    Returns
    -------
    (M,) ndarray
        g_n(0, r_2, ..., r_n) at unit number density.
    """
    m = r.shape[1]
    half = 0.5 * a
    # Breakpoints in u: r_j + u - a/2 crosses an integer (a per j inside [-a/2, a/2)).
    bps = -half + np.mod(-r, 1.0)[:, :, None] + np.arange(a)
    edges = np.sort(
        np.concatenate(
            [
                np.full((len(r), 1), -half),
                bps.reshape(len(r), -1),
                np.full((len(r), 1), half),
            ],
            axis=1,
        ),
        axis=1,
    )
    width = np.diff(edges, axis=1)
    mid = 0.5 * (edges[:, :-1] + edges[:, 1:])
    # Lattice sites k with |y - k| < a/2 are kmin, ..., kmin + a - 1.
    kmin = np.floor(r[:, None, :] + mid[:, :, None] - half) + 1
    count = np.zeros_like(width)
    for c in itertools.product(range(a), repeat=m):
        k = kmin + np.array(c)
        ok = np.all(k != 0, axis=-1)
        for i, j in itertools.combinations(range(m), 2):
            ok &= k[..., i] != k[..., j]
        count += ok
    # E_u carries density 1/a, and each neighbor's displacement density is 1/a.
    return np.sum(width * count, axis=1) / a ** (m + 1)


def bin_average(centers: NDArray, dx: float, a: int, sub: int) -> NDArray:
    """Average of the exact g_n over cubic bins of side dx (midpoint rule).

    Parameters
    ----------
    centers : (M, n-1) ndarray
        Bin centers [lattice constant].
    dx : float
        Bin width [lattice constant].
    a : int
        URL perturbation strength.
    sub : int
        Sub-samples per dimension and bin.

    Returns
    -------
    (M,) ndarray
        Bin-averaged g_n.
    """
    m = centers.shape[1]
    s = ((np.arange(sub) + 0.5) / sub - 0.5) * dx
    offsets = np.array(list(itertools.product(s, repeat=m)))
    pts = (centers[:, None, :] + offsets[None, :, :]).reshape(-1, m)
    return gn_url_exact(pts, a).reshape(len(centers), -1).mean(axis=1)


def load_gn(prefix: str, n: int) -> tuple[NDArray, NDArray, NDArray]:
    """Read the g_n text output of gn.out.

    Parameters
    ----------
    prefix : str
        Output prefix given to gn.out.
    n : int
        Order (2, 3, or 4).

    Returns
    -------
    centers : (M, n-1) ndarray
        Bin centers (x_21, ..., x_n1) [lattice constant].
    y, dy : (M,) ndarray
        g_n and its standard error (dy = 0 when absent).
    """
    if n in (2, 3):
        d = np.loadtxt(f"{prefix}_g{n}.txt", comments="#", ndmin=2)
        blocks = [d]
    elif n == 4:
        files = sorted(Path(prefix).parent.glob(Path(prefix).name + "_g4_x21_*.txt"))
        if not files:
            raise FileNotFoundError(f"no {prefix}_g4_x21_*.txt files")
        blocks = []
        for f in files:
            x21 = float(re.search(r"_g4_x21_([-0-9.]+)\.txt$", f.name).group(1))
            d = np.loadtxt(f, comments="#", ndmin=2)
            blocks.append(np.column_stack([np.full(len(d), x21), d]))
    else:
        raise ValueError(f"n must be 2, 3, or 4, got {n}")
    d = np.concatenate(blocks)
    centers = d[:, : n - 1]
    y = d[:, n - 1]
    dy = d[:, n] if d.shape[1] > n else np.zeros_like(y)
    return centers, y, dy


def main() -> int:
    """Compare gn.out output with the exact URL g_n (or with another run)."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("prefix", help="gn.out output prefix")
    ap.add_argument("--n", type=int, required=True, help="order of g_n (2, 3, 4)")
    ap.add_argument(
        "--a", type=int, default=1, help="URL perturbation strength (integer)"
    )
    ap.add_argument(
        "--sub", type=int, default=4, help="sub-samples per dimension per bin"
    )
    ap.add_argument(
        "--compare", help="second output prefix: report differences instead"
    )
    ap.add_argument(
        "--interior",
        action="store_true",
        help="with --compare: skip the outermost bins |i| = nh",
    )
    args = ap.parse_args()

    centers, y, dy = load_gn(args.prefix, args.n)
    if args.compare:
        c2, y2, dy2 = load_gn(args.compare, args.n)
        if c2.shape != centers.shape or not np.allclose(c2, centers):
            print("FAIL: bin grids differ")
            return 1
        sel = np.ones(len(y), dtype=bool)
        if args.interior:
            dx = float(np.min(np.diff(np.unique(centers[:, 0]))))
            nh = int(np.rint(centers[:, 0].max() / dx))
            sel = np.all(np.rint(np.abs(centers) / dx) < nh, axis=1)
        dmax = np.max(np.abs(y - y2)[sel])
        print(
            f"bins={sel.sum()}  max|y1-y2|={dmax:.3e}  "
            f"max|dy1-dy2|={np.max(np.abs(dy - dy2)[sel]):.3e}"
        )
        return 0 if dmax <= 1e-7 * max(1.0, np.max(np.abs(y))) else 1

    x21 = np.unique(centers[:, 0])  # x_21 >= 0 only
    dx = float(np.min(np.diff(x21))) if len(x21) > 1 else 1.0
    nh = int(np.rint(x21.max() / dx))
    ref = bin_average(centers, dx, args.a, args.sub)
    idx = np.rint(np.abs(centers) / dx).astype(int)
    edge = np.any(idx >= nh, axis=1)

    ok = True
    for name, sel in (("interior", ~edge), ("edge |i|=nh", edge)):
        if not np.any(sel):
            continue
        diff = y[sel] - ref[sel]
        has_err = dy[sel] > 0
        z = diff[has_err] / dy[sel][has_err]
        # A bin without samples (dy = 0) is expected where g_n is tiny: tuple counts there are
        # few and overdispersed (one n-particle cluster fills several bins), so flag only
        # empty bins where g_n is substantial.
        zero_err_bad = np.sum(~has_err & (np.abs(diff) > ZERO_ERR_TOL))
        line = (
            f"{name:12s} bins={sel.sum():6d}  mean(y-ref)={diff.mean():+.2e}  "
            f"max|y-ref|={np.max(np.abs(diff)):.3e}"
        )
        if z.size:
            line += (
                f"  rms z={np.sqrt(np.mean(z**2)):.3f}  max|z|={np.max(np.abs(z)):.2f}  "
                f"frac|z|>3={np.mean(np.abs(z) > 3):.4f}"
            )
        line += f"  dy=0: {np.sum(~has_err)} (|y-ref|>{ZERO_ERR_TOL}: {zero_err_bad})"
        print(line)
        ok &= zero_err_bad == 0
        if z.size:
            ok &= np.sqrt(np.mean(z**2)) < 1.3 and np.max(np.abs(z)) < 6
    print(f"dx={dx:g}  nh={nh}  a={args.a}  n={args.n}  ->  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
