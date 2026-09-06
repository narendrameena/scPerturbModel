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

Both reduce total sequencing. If the premise is right, thinning cells leaves the readout nearly unchanged until
the very lowest depths, while thinning conditions degrades it immediately. If
instead cells and conditions cost the same, the calculator is wrong and sec.46
should be withdrawn.

**The readout has to be a quantity Tahoe actually has.** A first version measured
the pooled interaction covariance, which on this atlas is -0.00016 -- indis-
tinguishable from zero, exactly as sec.31 found and as the design calculation
predicts. Thinning something that is already absent measures nothing. The readout
used instead is the effect sec.36 established IS present in Tahoe: MEK inhibitors
suppress the Pratilas ERK-output signature further in BRAF/RAS-driven lines than
in wild-type ones. It is large, biologically specified, and independently
validated, so its decay under each kind of thinning is interpretable.

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


MEK_OUTPUT = ["DUSP4", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1",
              "EPHA2", "SPRED1", "SPRED2", "CCND1", "FOSL1", "MYC"]
MEK_DRUGS = ["Cobimetinib", "Trametinib", "Binimetinib", "TAK-733"]
MAPK = ("BRAF", "KRAS", "NRAS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=12)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    C = pd.read_csv(PB / "conditions.csv")
    G = pd.read_csv(PB / "genes.csv")
    X = np.load(PB / "pseudobulk_counts.npz")["counts"]
    keep = (C.n_cells >= 200).to_numpy()
    C, X = C[keep].reset_index(drop=True), X[keep]
    sym = G.gene_symbol.astype(str).str.upper().to_numpy()
    gi = np.array([i for i, g in enumerate(sym) if g in set(MEK_OUTPUT)])
    md = pd.read_csv(TAB / "cell_line_metadata.csv")
    drv = md.groupby("Cell_ID_Cellosaur").Driver_Gene_Symbol.apply(
        lambda z: set(z.dropna().astype(str)))
    mapk = {k: bool(v & set(MAPK)) for k, v in drv.items()}
    print(f"{len(C)} pseudobulks; {len(gi)} ERK-output genes; "
          f"{sum(mapk.values())} MAPK-driven lines", flush=True)

    def prep(Xc):
        s_ = Xc.sum(1, keepdims=True)
        s_[s_ == 0] = 1.0
        return np.log1p(Xc / s_ * 1e4).astype(np.float32)

    def readout(Xl, cond, frac_cond=1.0, seed=0):
        """The sec.36 effect: how much further ERK output falls in MAPK lines."""
        r = np.random.default_rng(seed)
        ctl = {}
        for (ln, pl), g in cond[cond.drug.astype(str) == CTRL].groupby(
                ["cell_line_id", "plate"], observed=True):
            ctl[(ln, pl)] = Xl[g.index.to_numpy()].mean(0)
        sel = cond[cond.drug.isin(MEK_DRUGS)
                   & (cond.conc == cond.conc.max())]
        if frac_cond < 1.0 and len(sel):
            lines_u = sorted(set(sel.cell_line_id) & set(mapk))
            k = max(int(len(lines_u) * frac_cond), 4)
            keepl = set(r.choice(lines_u, k, replace=False))
            sel = sel[sel.cell_line_id.isin(keepl)]
        a, b = [], []
        for t in sel.itertuples():
            key = (t.cell_line_id, t.plate)
            if key not in ctl or t.cell_line_id not in mapk:
                continue
            d = Xl[t.Index][gi] - ctl[key][gi]
            (a if mapk[t.cell_line_id] else b).append(float(d.mean()))
        if len(a) < 5 or len(b) < 5:
            return np.nan
        return float(np.median(a) - np.median(b))

    rng = np.random.default_rng(0)
    rows = []
    print("\n1. THINNING CELLS (all conditions kept)", flush=True)
    for frac in (1.0, 0.5, 0.25, 0.1, 0.05, 0.02):
        vals = []
        for b_ in range(1 if frac >= 1.0 else args.n_boot // 3):
            Xt = X if frac >= 1.0 else rng.binomial(
                X.astype(np.int64), frac).astype(np.float32)
            v = readout(prep(Xt), C, 1.0, seed=b_)
            if np.isfinite(v):
                vals.append(v)
        if vals:
            rows.append({"axis": "cells", "frac": frac,
                         "effect": float(np.mean(vals)),
                         "sd": float(np.std(vals)) if len(vals) > 1 else 0.0})
            print(f"   {frac:5.0%} of cells: MAPK-vs-wildtype gap "
                  f"{np.mean(vals):+.4f}", flush=True)

    print("\n2. THINNING CONTEXTS (cells kept)", flush=True)
    Xf = prep(X)
    for frac in (1.0, 0.5, 0.25, 0.1):
        vals = []
        for b_ in range(1 if frac >= 1.0 else args.n_boot):
            v = readout(Xf, C, frac, seed=100 + b_)
            if np.isfinite(v):
                vals.append(v)
        if vals:
            rows.append({"axis": "contexts", "frac": frac,
                         "effect": float(np.mean(vals)),
                         "sd": float(np.std(vals)) if len(vals) > 1 else 0.0})
            print(f"   {frac:5.0%} of contexts: gap {np.mean(vals):+.4f} "
                  f"(sd {np.std(vals):.4f})", flush=True)

    T = pd.DataFrame(rows)
    T.to_csv(TAB / "cells_vs_replicates.csv", index=False)
    cel = T[T.axis == "cells"].set_index("frac").effect
    con = T[T.axis == "contexts"].set_index("frac").effect
    print("\n3. THE EXCHANGE RATE")
    base = cel.get(1.0, np.nan)
    if np.isfinite(base) and base != 0:
        for f in (0.1, 0.05, 0.02):
            if f in cel.index:
                print(f"   {f:.0%} of the cells retains "
                      f"{cel[f] / base:.0%} of the effect")
        for f in (0.25, 0.1):
            if f in con.index:
                print(f"   {f:.0%} of the contexts retains "
                      f"{con[f] / base:.0%} of the effect, with spread "
                      f"{T[(T.axis=='contexts')&(T.frac==f)].sd.iloc[0]:.4f}")
    print("   The premise of sec.46 predicts the first row to stay near 100% "
          "far down and\n   the second to degrade and become unstable. If both "
          "decay alike, cells and\n   contexts cost the same and sec.46 is "
          "wrong.")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    for axis, col, lab in (("cells", BLUE, "thin cells per condition"),
                           ("contexts", ORANGE, "thin contexts")):
        g = T[T.axis == axis].sort_values("frac")
        if len(g):
            ax[0].errorbar(g.frac, g.effect, yerr=g.sd, fmt="o-", color=col,
                           lw=2, ms=6, capsize=3, label=lab)
    ax[0].axhline(0, color="#444", lw=0.9)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("fraction of the data retained")
    ax[0].set_ylabel("MAPK-driven minus wild-type ERK suppression")
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title("a  A signal Tahoe demonstrably has (§36)", loc="left",
                    fontweight="bold", fontsize=9.5)

    if np.isfinite(base) and base != 0:
        for axis, col in (("cells", BLUE), ("contexts", ORANGE)):
            g = T[T.axis == axis].sort_values("frac")
            if len(g):
                ax[1].plot(g.frac, g.effect / base, "o-", color=col, lw=2,
                           ms=6, label=axis)
    ax[1].axhline(1.0, ls="--", color="#555", lw=1.3)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("fraction retained")
    ax[1].set_ylabel("effect relative to the full data")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].set_title("b  Which axis actually costs", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Testing the design calculation's premise on a signal that is "
                 "known to be present", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    d = save_figure(fig, "cells_vs_replicates", FIG, source_data={"sweep": T},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
