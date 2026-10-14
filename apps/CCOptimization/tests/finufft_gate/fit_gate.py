"""Fit the cost constants of ShiftedCCPotential (Type 2 force gate, Type 1 blocks) from calibrate_gate.out output.

Usage: python3 fit_gate.py CALIBRATION_FILE

Per dimension, the Type 2 time of one execution, type2 = Call + Particle * N + Mode * G [us], is the
tightest linear upper envelope of the measured times (non-negative coefficients, every case on or below
it, smallest summed relative prediction), rounded up to three significant digits; the Type 1 time of one block,
type1 = Call + Particle * n1 + Mode * G, is fitted the same way. Pair is the smallest
single-thread direct cost per particle-constraint pair among the long direct loops (short loops dodge
load and overstate the speed), rounded down. The script prints the four constant lines in the form of
CollectiveCoordinatePotential.h, the envelope's prediction/measurement ratios, and a margin check: with
the fitted constants, the cases the single-object gate accepts (direct >= MARGIN * type2, direct =
Pair * N * M / T) must be measured wins (direct loop wall time over T threads / Type 2 time).
"""
import itertools
import sys

import numpy as np

MARGIN = 1.3  # FinufftForceMargin in CollectiveCoordinatePotential.h
LONG_LOOP_PAIRS = 4e4  # direct loops with at least this many pairs define Pair
COMMENTS = ["Type 2 fixed cost per execution, by dimension.",
            "Type 2 cost per particle (all d transforms) at kernel width 13.",
            "Type 2 cost per mode-grid point.",
            "Direct Force() cost per particle-constraint pair.",
            "Type 1 fixed cost per execution, by dimension.",
            "Type 1 cost per particle at kernel width 13.",
            "Type 1 cost per mode-grid point."]


def load(path):
    """Read calibrate_gate.out lines.

    Parameters
    ----------
    path : str
        File with lines ``skew chi dim N M G type2_us pair_us directT1_us directT2_us n1 type1_us``; lines
        starting with ``#`` hold comments and the two team sizes T1, T2 in the header.

    Returns
    -------
    rows : numpy.ndarray, shape (cases, 12)
        One row per case, columns as above; times in microseconds.
    teams : tuple of int
        Team sizes of the two direct wall-time columns.
    """
    teams = (4, 8)
    with open(path) as handle:
        for line in handle:
            if line.startswith("# skew"):
                teams = tuple(int(word[len("direct"):-len("_us")]) for word in line.split() if word.startswith("direct"))
    rows = np.loadtxt(path, comments="#", ndmin=2)
    return rows, teams


def upper_envelope(n, g, t):
    """Tightest linear upper envelope t <= a + b n + c g with a, b, c >= 0.

    Solves the linear program min sum_i (a + b n_i + c g_i) / t_i subject to the envelope constraints by
    enumerating the vertices of the feasible set (three active constraints among the cases and a, b, c = 0).

    Parameters
    ----------
    n, g, t : numpy.ndarray, shape (cases,)
        Particle numbers, mode-grid sizes and measured times [us].

    Returns
    -------
    numpy.ndarray, shape (3,)
        Coefficients (a [us], b [us per particle], c [us per mode]).
    """
    scale = np.array([1.0, n.max(), g.max()])
    a = np.column_stack([np.ones_like(n), n, g]) / scale  # columns of order one
    rows = np.vstack([a, np.eye(3)])
    rhs = np.concatenate([t, np.zeros(3)])
    trio = np.array(list(itertools.combinations(range(len(rows)), 3)))
    m, b = rows[trio], rhs[trio]
    regular = np.abs(np.linalg.det(m)) > 1e-12
    x = np.linalg.solve(m[regular], b[regular][..., None])[..., 0]
    feasible = (x >= -1e-12).all(axis=1) & (a @ x.T >= t[:, None] * (1 - 1e-9)).all(axis=0)
    cost = (a / t[:, None]).sum(axis=0)
    best = x[feasible][np.argmin(x[feasible] @ cost)]
    return np.maximum(best, 0.0) / scale


def round_sig(x, up):
    """Round to three significant digits, up or down.

    Parameters
    ----------
    x : float
        Non-negative value.
    up : bool
        Round up (True) or down (False).

    Returns
    -------
    float
        The rounded value.
    """
    if x <= 0 or not np.isfinite(x):
        return x
    step = 10.0 ** (np.floor(np.log10(x)) - 2)
    return float((np.ceil if up else np.floor)(x / step - (1e-9 if up else -1e-9)) * step)


def main():
    """Fit, print the constants and the checks; return 0, or 1 on bad usage."""
    if len(sys.argv) != 2:
        print(__doc__)
        return 1
    rows, teams = load(sys.argv[1])
    skew, chi, dim, n, m, g, t2, pair, wall1, wall2, n1, t1 = rows.T
    consts = np.full((7, 3), np.nan)
    print("# envelope prediction/measurement per dimension (min must be 1.00):")
    for d in (1, 2, 3):
        sel = dim == d
        if sel.sum() < 3:
            print(f"#   {d}D: fewer than 3 cases, constants left as nan")
            continue
        coef = upper_envelope(n[sel], g[sel], t2[sel])
        consts[:3, d - 1] = [round_sig(c, up=True) for c in coef]
        ratio = (coef[0] + coef[1] * n[sel] + coef[2] * g[sel]) / t2[sel]
        long_loops = sel & (n * m >= LONG_LOOP_PAIRS)
        rates = pair[long_loops] if long_loops.any() else pair[sel]
        consts[3, d - 1] = round_sig(rates.min(), up=False)
        coef1 = upper_envelope(n1[sel], g[sel], t1[sel])
        consts[4:, d - 1] = [round_sig(c, up=True) for c in coef1]
        ratio1 = (coef1[0] + coef1[1] * n1[sel] + coef1[2] * g[sel]) / t1[sel]
        print(f"#   {d}D: {sel.sum()} cases, Type 2 min {ratio.min():.2f} median {np.median(ratio):.2f} max {ratio.max():.2f}, "
              f"Type 1 min {ratio1.min():.2f} median {np.median(ratio1):.2f} max {ratio1.max():.2f}")
    names = ["FinufftCostCall", "FinufftCostParticle", "FinufftCostMode", "FinufftCostPair",
             "FinufftRhoCostCall", "FinufftRhoCostParticle", "FinufftRhoCostMode"]
    for name, values, comment in zip(names, consts, COMMENTS):
        text = ", ".join(f"{v:.3g}" for v in values)
        print(f"\tstatic constexpr double {name}[3]={{{text}}};\t//< {comment}")

    print(f"# margin check with the fitted constants (single object, eps 1e-12, margin {MARGIN}):")
    idx = (dim - 1).astype(int)
    model = consts[0, idx] + consts[1, idx] * n + consts[2, idx] * g
    fitted = np.isfinite(model)
    for team, wall in zip(teams, (wall1, wall2)):
        accepted = fitted & (consts[3, idx] * n * m / team >= MARGIN * model)
        speedup = wall / t2
        missed = fitted & (~accepted) & (speedup >= MARGIN)
        line = f"#   T={team}: {accepted.sum()} of {fitted.sum()} fitted cases accepted"
        if accepted.any():
            line += f", measured speedup min {speedup[accepted].min():.2f} median {np.median(speedup[accepted]):.2f}"
            if speedup[accepted].min() < 1.0:
                line += "  WARNING: an accepted case is slower than the direct loop"
        print(line + f"; {missed.sum()} declined cases measured a speedup >= {MARGIN}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
