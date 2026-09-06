#!/usr/bin/env python3
"""Cells or replicates? The exchange rate, measured rather than assumed.

RESULTS.md sec.46 built a design calculation whose central premise is that the
precision of an interaction estimate is set by the number of replicate PAIRS and
not by the number of cells. That premise is the reason the calculator tells
atlas builders to spend differently, so it should be tested on real data rather
than asserted from the algebra — and it is the kind of claim that can fail.

The honest version of the premise is not "cells do not matter". Cells enter
through the per-condition noise: halve the cells and each condition's measurement
gets noisier, which enters the estimator's variance through the signal-to-noise
term. What the algebra says is that this dependence SATURATES — past the point
where a condition's measurement is dominated by biological rather than sampling
variation, more cells buy almost nothing — while each added replicate creates
new pairs and buys precision linearly.

So the experiment is a two-way subsampling of Tahoe:

  * cells per condition, thinned from 100% down to 5%, replicates held fixed;
  * conditions, thinned from 100% down to 5%, cells held fixed.

Both reduce total sequencing. If the premise is right, thinning cells leaves the
interaction estimate and its standard error nearly unchanged until the very
lowest depths, while thinning conditions degrades both immediately and
predictably. If instead cells and conditions cost the same, the calculator is
wrong and sec.46 should be withdrawn.

The measured exchange rate then answers the practical question directly: for a
fixed budget of sequenced cells, how should they be spread?

Outputs: results/tables/cells_vs_replicates.csv
         figure bundle results/figures/00_manuscript/cells_vs_replicates/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PB = ROOT / "data" / "processed" / "pseudobulk_full"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
CTRL = "DMSO_TF"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=30)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    C = pd.read_csv(PB / "conditions.csv")
    X = np.load(PB / "pseudobulk_counts.npz")["counts"]
    keep = (C.n_cells >= 200).to_numpy()
    C, X = C[keep].reset_index(drop=True), X[keep]
    print(f"{len(C)} pseudobulks", flush=True)

    # Thinning cells is simulated by binomial downsampling of the counts, which
    # is what sequencing fewer cells does to a pseudobulk: the same expected
    # profile with sampling noise scaled by 1/sqrt(depth). Re-normalising after
    # thinning keeps the scale comparable.
    def prep(Xc):
        s = Xc.sum(1, keepdims=True)
        s[s == 0] = 1.0
        return np.log1p(Xc / s * 1e4).astype(np.float32)

    rng = np.random.default_rng(0)
    var = np.asarray(X.sum(0)).ravel()
    resp = np.argsort(-var)[:2000]

    def interaction(Xl, cond, frac_cond=1.0, seed=0):
        """Estimate the interaction share on a (possibly thinned) dataset."""
        r = np.random.default_rng(seed)
        ctl = {}
        for (ln, pl), g in cond[cond.drug.astype(str) == CTRL].groupby(
                ["cell_line_id", "plate"], observed=True):
            ctl[(ln, pl)] = Xl[g.index.to_numpy()].mean(0)
        trt = cond[cond.drug.astype(str) != CTRL]
        rows = [(t.cell_line_id, t.drug, t.conc, t.plate, i)
                for i, t in zip(trt.index.to_numpy(), trt.itertuples())
                if (t.cell_line_id, t.plate) in ctl]
        K = pd.DataFrame(rows, columns=["line", "drug", "conc", "plate", "i"])
        if frac_cond < 1.0:
            keys = K.groupby(["drug", "conc"], observed=True).ngroup()
            uq = np.unique(keys)
            pick = set(r.choice(uq, max(int(len(uq) * frac_cond), 2),
                                replace=False))
            K = K[keys.isin(pick)].reset_index(drop=True)
        if not len(K):
            return np.nan
        D = np.stack([Xl[t.i] - ctl[(t.line, t.plate)] for t in K.itertuples()])
        # drug main effect out, leave-one-line-out per (drug, dose)
        res = np.zeros_like(D)
        for (dr, cc), g in K.groupby(["drug", "conc"], observed=True):
            ii = g.index.to_numpy()
            ln = K.line.to_numpy()[ii]
            if len(np.unique(ln)) < 2:
                continue
            tot = D[ii].sum(0)
            cs = {c: D[ii[ln == c]].sum(0) for c in np.unique(ln)}
            cn = {c: int((ln == c).sum()) for c in np.unique(ln)}
            for i, c in zip(ii, ln):
                n_out = len(ii) - cn[c]
                if n_out >= 1:
                    res[i] = D[i] - (tot - cs[c]) / n_out
        # line general response out, within plate
        dv = K.drug.to_numpy()
        for (ln, pl), g in K.groupby(["line", "plate"], observed=True):
            ii = g.index.to_numpy()
            by = {}
            for i in ii:
                by.setdefault(dv[i], []).append(res[i])
            by = {d: np.mean(v, axis=0) for d, v in by.items()}
            if len(by) < 3:
                continue
            tot_a, n_a = np.sum(list(by.values()), axis=0), len(by)
            for i in ii:
                res[i] = res[i] - (tot_a - by[dv[i]]) / (n_a - 1)
        # cross-plate covariance at matched (line, drug, dose)
        num, den, n = 0.0, 0.0, 0
        for (ln, dr, cc), g in K.groupby(["line", "drug", "conc"],
                                         observed=True):
            ii = g.i.to_numpy()
            pl = K.plate.to_numpy()[g.index.to_numpy()]
            if len(set(pl)) < 2:
                continue
            a = g.index.to_numpy()[pl == pl[0]]
            b = g.index.to_numpy()[pl != pl[0]]
            num += float(np.mean(res[a].mean(0) * res[b].mean(0))); n += 1
        if n < 20:
            return np.nan
        return num / n

    base_cond = C.copy()
    print("\n1. THINNING CELLS (replicates and conditions held fixed)",
          flush=True)
    rows = []
    for frac in (1.0, 0.5, 0.25, 0.1, 0.05):
        vals = []
        for b in range(args.n_boot // 6 if frac < 1 else 1):
            if frac >= 1.0:
                Xt = X
            else:
                Xt = rng.binomial(X.astype(np.int64),
                                  frac).astype(np.float32)
            v = interaction(prep(Xt)[:, resp], base_cond, 1.0, seed=b)
            if np.isfinite(v):
                vals.append(v)
        if vals:
            rows.append({"axis": "cells", "frac": frac,
                         "estimate": float(np.mean(vals)),
                         "sd": float(np.std(vals)) if len(vals) > 1 else 0.0,
                         "n": len(vals)})
            print(f"   {frac:5.0%} of cells: interaction "
                  f"{np.mean(vals):.5f}", flush=True)

    print("\n2. THINNING CONDITIONS (cells held fixed)", flush=True)
    Xf = prep(X)[:, resp]
    for frac in (1.0, 0.5, 0.25, 0.1, 0.05):
        vals = []
        for b in range(1 if frac >= 1.0 else args.n_boot // 3):
            v = interaction(Xf, base_cond, frac, seed=100 + b)
            if np.isfinite(v):
                vals.append(v)
        if vals:
            rows.append({"axis": "conditions", "frac": frac,
                         "estimate": float(np.mean(vals)),
                         "sd": float(np.std(vals)) if len(vals) > 1 else 0.0,
                         "n": len(vals)})
            print(f"   {frac:5.0%} of conditions: interaction "
                  f"{np.mean(vals):.5f}  (sd {np.std(vals):.5f} over "
                  f"{len(vals)} draws)", flush=True)

    T = pd.DataFrame(rows)
    T.to_csv(TAB / "cells_vs_replicates.csv", index=False)

    cel = T[T.axis == "cells"].set_index("frac").estimate
    con = T[T.axis == "conditions"].set_index("frac").estimate
    print("\n3. THE EXCHANGE RATE")
    if 1.0 in cel.index and 0.1 in cel.index:
        print(f"   cutting cells to 10% changes the estimate by "
              f"{abs(cel[0.1] - cel[1.0]) / abs(cel[1.0]):.0%}")
    if 1.0 in con.index and 0.1 in con.index:
        print(f"   cutting conditions to 10% changes it by "
              f"{abs(con[0.1] - con[1.0]) / abs(con[1.0]):.0%}, and its "
              f"spread across draws by "
              f"{T[(T.axis=='conditions')&(T.frac==0.1)].sd.iloc[0]:.5f}")
    print("   If the first is small and the second large, the premise of "
          "sec.46 holds:\n   sequencing depth per condition is not the binding "
          "constraint, and the same\n   cells spread over more conditions and "
          "replicates buy more information.")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    for axis, col, lab in (("cells", BLUE, "thin the cells per condition"),
                           ("conditions", ORANGE, "thin the conditions")):
        g = T[T.axis == axis].sort_values("frac")
        if len(g):
            ax[0].errorbar(g.frac, g.estimate, yerr=g.sd, fmt="o-", color=col,
                           lw=2, ms=6, capsize=3, label=lab)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("fraction of the data retained")
    ax[0].set_ylabel("estimated interaction (cross-plate covariance)")
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title("a  Two ways to spend less sequencing", loc="left",
                    fontweight="bold", fontsize=9.5)

    base_c = cel.get(1.0, np.nan)
    for axis, col in (("cells", BLUE), ("conditions", ORANGE)):
        g = T[T.axis == axis].sort_values("frac")
        if len(g) and np.isfinite(base_c) and base_c != 0:
            ax[1].plot(g.frac, (g.estimate / base_c), "o-", color=col, lw=2,
                       ms=6, label=axis)
    ax[1].axhline(1.0, ls="--", color="#555", lw=1.3)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("fraction retained")
    ax[1].set_ylabel("estimate relative to the full data")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].text(0.5, 0.1, "flat under cell thinning = depth is not\nthe binding "
               "constraint", transform=ax[1].transAxes, ha="center",
               fontsize=7, color="#444")
    ax[1].set_title("b  Which axis actually costs", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Testing the design calculation's premise: cells enter only "
                 "through per-condition noise", fontsize=10.5, x=0.005,
                 ha="left", fontweight="bold")
    d = save_figure(fig, "cells_vs_replicates", FIG, source_data={"sweep": T},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
