#!/usr/bin/env python3
"""Are contexts and replicates interchangeable at equal pair count? PRISM says.

This is the test the paper's title makes and the paper never ran. The design
calculation asserts that three design numbers enter the precision of an
interaction estimate through exactly one combination:

    n_pairs = n_ctx x n_pert x n_rep(n_rep-1)/2      and SE ~ n_pairs^(-1/2)

Everything prescriptive in the paper follows from that: that a replicate and a
context are worth the same at equal pair count, that cells are worth nothing, and
that a budget should be spread rather than deepened. RESULTS.md sec.46 validated
the *verdicts* the formula produces on five atlases; the 2026-09-10 audit showed
those verdicts are reproduced by any constant floor in a 21-fold window, so they
do not test this. Nothing else does either.

**Design.** PRISM's secondary screen is a 737 x 1,488 x 3 cube (lines x compounds
x replicate detection plates) with a *scalar* readout, so `n_feat` drops out and
the design terms are isolated. Sub-cubes are drawn at a grid of
(n_ctx, n_pert, n_rep), the interaction share is estimated on each, and the
standard error is measured as the spread across draws. Then:

    log SE ~ a.log(n_ctx) + b.log(n_pert) + c.log(n_rep(n_rep-1)/2)

The calculation asserts **a = b = c = -1/2**. Three things can fail, each
informative:

  * `c != -1/2`  -- replicates do not buy what the formula says they buy.
  * `a != c` or `b != c` -- contexts and replicates are NOT interchangeable at
    equal pair count. This is the one most likely to fail, because a second
    replicate of the same context shares that context's biology while a new
    context does not, and it is the claim the title makes.
  * the level `SE_measured / SE_predicted` not constant across the grid --
    the formula has the shape right but no calibration.

**What this does not do.** It never resamples pairs. The withdrawn sec.49 test
bootstrapped pair-products i.i.d., which forces SE ~ n_pairs^(-1/2) as an
identity and so could not have measured this. Here the resampling unit is the
line and the compound, and the pair structure is rebuilt inside each draw.

**Control.** The whole grid is rerun with compound labels permuted within line,
destroying the interaction while preserving both marginals. The SE should obey
the same design law while the estimate collapses; if the law changes, the
machinery is tracking signal rather than variance.

Outputs: results/tables/exchange_rate.csv
         figure bundle results/figures/00_manuscript/exchange_rate/
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"

N_CTX = [5, 10, 20, 40, 80, 160]
N_PERT = [10, 25, 50, 100, 200, 400]
N_REP = [2, 3]
MIN_FINITE = 0.30      # a drawn sub-cube needs this fraction of finite cells


def build_cube(R, K):
    """(compound,dose) x rep -> vector over lines, as a dense 3-D array.

    Rows are conditions with all three replicates present, so n_rep can be varied
    by slicing rather than by re-selecting conditions -- otherwise changing n_rep
    would also change which conditions are in play and confound the comparison.
    """
    reps = sorted(K.rep.unique())
    idx = {}
    for i, (c, d, r) in enumerate(zip(K.compound, K.dose, K.rep)):
        idx.setdefault((c, d), {})[r] = i
    keep = [(cd, m) for cd, m in idx.items() if all(r in m for r in reps)]
    cube = np.full((len(keep), len(reps), R.shape[0]), np.nan, dtype=np.float32)
    meta = []
    for j, (cd, m) in enumerate(keep):
        for k, r in enumerate(reps):
            cube[j, k] = R[:, m[r]]
        meta.append(cd)
    return cube, pd.DataFrame(meta, columns=["compound", "dose"]), reps


def interaction_share(sub):
    """Interaction share of a (conditions, reps, lines) sub-cube.

    Both main effects out of fold: the condition effect from other lines, the
    line effect from other conditions. The interaction is the covariance between
    independent replicates, so noise contributes zero in expectation.
    """
    C, Rn, L = sub.shape
    if C < 3 or L < 3 or Rn < 2:
        return np.nan
    m = np.nanmean(sub, axis=1)                      # conditions x lines
    ok = np.isfinite(m)
    if ok.mean() < MIN_FINITE:
        return np.nan
    # condition effect, leave-one-line-out
    tot = np.nansum(m, axis=1, keepdims=True)
    cnt = ok.sum(axis=1, keepdims=True)
    beta = (tot - np.nan_to_num(m)) / np.maximum(cnt - ok, 1)
    Rz = m - beta
    # line effect, leave-one-condition-out
    tot2 = np.nansum(Rz, axis=0, keepdims=True)
    cnt2 = np.isfinite(Rz).sum(axis=0, keepdims=True)
    alpha = (tot2 - np.nan_to_num(Rz)) / np.maximum(cnt2 - np.isfinite(Rz), 1)
    resid = Rz - alpha
    # residual per replicate, using the same main effects
    per = sub - beta[:, None, :] - alpha[:, None, :]
    cov, n = 0.0, 0
    for a in range(Rn):
        for b in range(a + 1, Rn):
            v = per[:, a, :] * per[:, b, :]
            f = np.isfinite(v)
            if f.sum():
                cov += float(np.nansum(v)); n += int(f.sum())
    if not n:
        return np.nan
    inter = max(cov / n, 0.0)
    shared = float(np.nanmean(resid ** 2))
    tot_ = inter + shared
    return inter / tot_ if tot_ > 0 else np.nan


def run_grid(cube, n_lines, rng, n_draws, permute=False):
    rows = []
    C = cube.shape[0]
    for nc in N_CTX:
        if nc > n_lines:
            continue
        for npt in N_PERT:
            if npt > C:
                continue
            for nr in N_REP:
                vals = []
                for _ in range(n_draws):
                    li = rng.choice(n_lines, nc, replace=False)
                    ci = rng.choice(C, npt, replace=False)
                    sub = cube[np.ix_(ci, np.arange(nr), li)]
                    if permute:
                        # Break the LINE x CONDITION interaction while keeping
                        # both main effects intact: permute the line axis
                        # independently within each replicate. Replicate a of a
                        # condition is then paired with a different line's
                        # replicate b, so the true interaction cannot survive,
                        # but every condition still has its shared effect and
                        # every line its general sensitivity.
                        #
                        # NOT permuting conditions within a line, which was the
                        # first attempt: that destroys the SHARED condition
                        # effect, leaving it inside the residual where both
                        # replicates still carry it, and drives the estimate UP
                        # (0.47 against 0.38 for real data) instead of to zero.
                        sub = sub.copy()
                        for k in range(sub.shape[1]):
                            sub[:, k, :] = sub[:, k, rng.permutation(sub.shape[2])]
                    v = interaction_share(sub)
                    if np.isfinite(v):
                        vals.append(v)
                if len(vals) < 8:
                    continue
                rows.append({"n_ctx": nc, "n_pert": npt, "n_rep": nr,
                             "pairs_per_cond": nr * (nr - 1) // 2,
                             "n_pairs": nc * npt * nr * (nr - 1) // 2,
                             "mean_share": float(np.mean(vals)),
                             "se": float(np.std(vals, ddof=1)),
                             "n_draws": len(vals),
                             "permuted": permute})
    return rows


def fit(D, label):
    """log SE on the three design terms; the calculation says all three are -0.5."""
    d = D[(D.se > 0) & np.isfinite(D.se)]
    X = np.column_stack([np.log(d.n_ctx), np.log(d.n_pert),
                         np.log(d.pairs_per_cond), np.ones(len(d))])
    y = np.log(d.se.to_numpy())
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    dof = len(d) - X.shape[1]
    s2 = float(resid @ resid) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    names = ["a (contexts)", "b (perturbations)", "c (replicate pairs)"]
    print(f"\n  {label}: log SE ~ a.log(n_ctx) + b.log(n_pert) "
          f"+ c.log(pairs/cond)")
    print(f"    n = {len(d)} grid cells, residual sd {np.sqrt(s2):.3f}")
    out = {}
    for i, nm in enumerate(names):
        lo, hi = coef[i] - 1.96 * se[i], coef[i] + 1.96 * se[i]
        hit = "contains -0.5" if lo <= -0.5 <= hi else "EXCLUDES -0.5"
        print(f"    {nm:22s} {coef[i]:+.3f} [{lo:+.3f}, {hi:+.3f}]  {hit}")
        out[nm] = (coef[i], lo, hi)
    # are contexts and replicates interchangeable? a == c
    diff = coef[0] - coef[2]
    sd = float(np.sqrt(cov[0, 0] + cov[2, 2] - 2 * cov[0, 2]))
    z = diff / sd if sd > 0 else np.nan
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    print(f"    a - c = {diff:+.3f} +/- {sd:.3f}   z = {z:+.2f}, P = {p:.4f}"
          f"   {'INTERCHANGEABLE' if p > 0.05 else 'NOT interchangeable'}")
    out["a_minus_c"] = (diff, sd, p)
    return out, coef


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-draws", type=int, default=30)
    a = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    from perturbmodel.celldrug import load_prism
    print("loading PRISM ...", flush=True)
    R, K, lines, _ = load_prism()
    cube, meta, reps = build_cube(R, K)
    print(f"  cube: {cube.shape[0]:,} conditions x {len(reps)} reps x "
          f"{cube.shape[2]} lines", flush=True)

    rng = np.random.default_rng(0)
    print(f"\nsweeping the design grid ({a.n_draws} draws per cell) ...",
          flush=True)
    rows = run_grid(cube, len(lines), rng, a.n_draws, permute=False)
    print(f"  {len(rows)} usable grid cells", flush=True)
    print("control: same grid with the interaction destroyed ...", flush=True)
    rows += run_grid(cube, len(lines), np.random.default_rng(1),
                     max(a.n_draws // 2, 10), permute=True)
    D = pd.DataFrame(rows)
    D.to_csv(TAB / "exchange_rate.csv", index=False)

    real = D[~D.permuted]
    perm = D[D.permuted]
    print(f"\n  mean interaction share: real {real.mean_share.mean():.4f}, "
          f"permuted {perm.mean_share.mean():.4f}"
          if len(perm) else "")

    out, coef = fit(real, "REAL DATA (raw SE)")
    # The CV fit MUST be reported per-term, not only on the aggregate. Reporting
    # it on the aggregate alone is how the withdrawn interchangeability claim
    # survived to publication: a - c is P = 0.016 on raw SE and P = 0.70 on CV,
    # and the estimand moves 7.8% with n_rep -- larger than the context drift the
    # aggregate check was controlling for.
    real_cv = real.assign(se=real.se / real.mean_share)
    fit(real_cv, "REAL DATA (coefficient of variation -- the drift control)")
    if len(perm) > 6:
        # Per-term, not just the mean-share collapse. This docstring's own
        # criterion is that the LAW should not change; only reporting the
        # collapse hid that it does.
        fit(perm, "PERMUTED CONTROL (interaction destroyed)")

    # level: is SE/predicted a single constant across the grid?
    pred = 1.0 / np.sqrt(real.n_pairs)
    ratio = real.se / pred
    print(f"\n  level check: SE / (1/sqrt(pairs)) spans "
          f"{ratio.min():.3f} to {ratio.max():.3f} "
          f"({ratio.max()/max(ratio.min(),1e-9):.1f}x), median {ratio.median():.3f}")

    print("\n  VERDICT")
    ca, la, ha = out["c (replicate pairs)"]
    print(f"    replicate exponent -0.5?      "
          f"{'YES' if la <= -0.5 <= ha else 'NO'}")
    _, _, p_int = out["a_minus_c"]
    print(f"    contexts == replicates?       "
          f"{'YES' if p_int > 0.05 else 'NO'}  (P = {p_int:.4f})")
    print(f"    level constant within 3x?     "
          f"{'YES' if ratio.max()/max(ratio.min(),1e-9) < 3 else 'NO'}")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True)
    for nr, col in zip(N_REP, (BLUE, ORANGE)):
        g = real[real.n_rep == nr]
        ax[0].scatter(g.n_pairs, g.se, s=26, color=col, alpha=0.85,
                      edgecolor="white", lw=0.5, label=f"{nr} replicates")
    xs = np.logspace(np.log10(real.n_pairs.min()), np.log10(real.n_pairs.max()), 50)
    k = float(np.median(real.se * np.sqrt(real.n_pairs)))
    ax[0].plot(xs, k / np.sqrt(xs), color="#111", lw=1.6, ls="--",
               label="predicted slope −0.5")
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("replicate pairs in the design")
    ax[0].set_ylabel("measured SE of the interaction share")
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title("a  Does precision follow the pair count?", loc="left",
                    fontweight="bold", fontsize=9.5)

    names = ["contexts", "perturbations", "replicate pairs"]
    vals = [out[f"{n} " if False else k][0] for k, n in
            zip(["a (contexts)", "b (perturbations)", "c (replicate pairs)"],
                names)]
    los = [out[k][1] for k in ["a (contexts)", "b (perturbations)",
                               "c (replicate pairs)"]]
    his = [out[k][2] for k in ["a (contexts)", "b (perturbations)",
                               "c (replicate pairs)"]]
    y = np.arange(3)[::-1]
    ax[1].errorbar(vals, y, xerr=[np.array(vals) - np.array(los),
                                  np.array(his) - np.array(vals)],
                   fmt="o", color=VIOLET, ms=8, capsize=4, lw=2)
    ax[1].axvline(-0.5, color=ORANGE, lw=2, ls="--", label="predicted −0.5")
    ax[1].set_yticks(y, names, fontsize=9)
    ax[1].set_xlabel("exponent on log SE")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].set_title("b  Are the three terms interchangeable?", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Testing the exchange rate the design calculation asserts, on "
                 "PRISM", fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "exchange_rate", FIG, source_data={"grid": D},
                    script=__file__)
    print(f"\nfigure bundle -> {d}")


if __name__ == "__main__":
    main()
