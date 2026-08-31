#!/usr/bin/env python3
"""Is our context-dependence index the right metric? Three alternatives compared.

Our CDI is a variance ratio: reproducible interaction / (conserved +
interaction). The published literature mostly does something different — it
correlates a compound's response PROFILE between contexts (Niepel 2017;
Subramanian/CMap 2017 report the fraction of compounds with similar signatures
panel-wide; TRADE correlates perturbation profiles between cell types). Before
resting a claim on our metric we should know whether it agrees with theirs.

Three metrics per compound, all computed on the same conditions and genes:

  cdi          ours. reproducible interaction variance / total reproducible
               variance. Pools all contexts; separates signal from noise using
               replicates; but the numerator is a covariance that can go
               negative and is then clamped, which loses resolution at the
               conserved end.

  r_cross      the published approach: mean correlation of the compound's
               response profile between different CONTEXTS. Low = context
               specific. Simple and familiar, but confounds context-specificity
               with the compound's own reproducibility — a weak or noisy
               compound scores low for the wrong reason.

  transfer     correction for attenuation: r_cross / r_replicate, where
               r_replicate is the same correlation between independent
               REPLICATES of the same context x compound. This is the
               reliability ceiling, so the ratio asks how much of what is
               reproducible actually transfers. 1 = everything reproducible
               transfers; 0 = nothing does.

`transfer` is arguably the most defensible: it uses the familiar correlation
scale of the published work while removing the confound that makes raw
correlation misleading. If all three rank mechanisms alike, the choice does not
matter and we can report the familiar one; if they differ, that is worth
knowing before publication.

Outputs: results/tables/context_metric_comparison.csv
         figure bundle results/figures/10_drugs/context_metric_comparison/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.evaluation.delta_eval import (build_deltas, load_pseudobulk,
                                                responsive_genes)
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "10_drugs"
TAB = ROOT / "results" / "tables"
N_RESP = 2000
MIN_PAIRS = 10
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def cor(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    args = ap.parse_args()

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond, keep_plate=True)
    resp = responsive_genes(DELTA, np.ones(len(G), bool))
    D = DELTA[:, resp].astype(np.float32)
    del DELTA
    K = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                      "drug": G.drug, "conc": G.conc, "plate": G.plate})
    print(f"{len(K)} conditions, {len(resp)} genes", flush=True)

    recs = []
    for drug, gd in K.groupby("drug", observed=True):
        if gd.line.nunique() < 10:
            continue
        # --- ours: variance decomposition ---
        # leave-one-context-out shared response: an in-sample mean forces the
        # residuals of n conditions to sum to zero, biasing their covariance to
        # -sigma^2/n even with no interaction present
        # leave-one-CONTEXT-out: the covariance below is between two plates
        # of the SAME line, so the subtracted mean must exclude that line
        # entirely. Leaving out only the single condition leaves its sibling
        # replicate in the mean and biases the covariance to about
        # -2 sigma^2/(n-1); the in-sample mean biases it to -sigma^2/n.
        conserved, resid = 0.0, {}
        for _, gc in gd.groupby("conc", observed=True):
            ii = gc.i.to_numpy(); ln = gc.line.to_numpy()
            if len(np.unique(ln)) < 2:
                continue
            tot = D[ii].sum(0)
            csum = {c: D[ii[ln == c]].sum(0) for c in np.unique(ln)}
            ccnt = {c: int((ln == c).sum()) for c in np.unique(ln)}
            for i, c in zip(ii, ln):
                n_out = len(ii) - ccnt[c]
                if n_out < 1:
                    continue
                loo = (tot - csum[c]) / n_out
                conserved += float(np.mean(loo ** 2))
                resid[i] = D[i] - loo
        conserved /= max(len(resid), 1)
        cov, npair = 0.0, 0
        for _, gc in gd.groupby("line", observed=True):
            v = gc.i.to_numpy(); pl = gc.plate.to_numpy()
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if pl[a] != pl[b]:
                        cov += float(np.mean(resid[v[a]] * resid[v[b]])); npair += 1
        if npair < MIN_PAIRS:
            continue
        context = max(cov / npair, 0.0)
        den = conserved + context
        cdi = context / den if den > 0 else np.nan

        # --- published: profile correlation between CONTEXTS ---
        rc = []
        for _, gc in gd.groupby("conc", observed=True):
            v = gc.i.to_numpy(); ln = gc.line.to_numpy(); pl = gc.plate.to_numpy()
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if ln[a] != ln[b] and pl[a] != pl[b]:
                        rc.append(cor(D[v[a]], D[v[b]]))
        # --- reliability: correlation between REPLICATES of the same pair ---
        rr = []
        for _, gc in gd.groupby(["line", "conc"], observed=True):
            v = gc.i.to_numpy(); pl = gc.plate.to_numpy()
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if pl[a] != pl[b]:
                        rr.append(cor(D[v[a]], D[v[b]]))
        if len(rc) < MIN_PAIRS or len(rr) < 3:
            continue
        r_cross = float(np.nanmedian(rc))
        r_rep = float(np.nanmedian(rr))
        transfer = r_cross / r_rep if r_rep > 0.05 else np.nan
        recs.append({"drug": drug, "n_lines": gd.line.nunique(),
                     "cdi": cdi, "r_cross": r_cross, "r_replicate": r_rep,
                     "transfer": transfer, "n_pairs": npair})

    res = pd.DataFrame(recs)
    dm = pd.read_parquet(ROOT / "data/metadata/metadata/drug_metadata.parquet")
    res = res.merge(dm[["drug", "moa-fine"]], on="drug", how="left")
    res.to_csv(TAB / "context_metric_comparison.csv", index=False)
    print(f"\n{len(res)} compounds with all three metrics")
    print(f"  median cdi {res.cdi.median():.3f}, r_cross {res.r_cross.median():.3f}, "
          f"r_replicate {res.r_replicate.median():.3f}, "
          f"transfer {res.transfer.median():.3f}")

    v = res.dropna(subset=["cdi", "r_cross", "transfer"])
    print("\npairwise agreement between metrics (Spearman, per compound):")
    for a, b in (("cdi", "r_cross"), ("cdi", "transfer"), ("r_cross", "transfer")):
        rho = stats.spearmanr(v[a], v[b])
        print(f"  {a:9s} vs {b:9s} rho={rho.statistic:+.3f}  p={rho.pvalue:.1e}")
    rho_eff = stats.spearmanr(v.r_cross, v.r_replicate)
    print(f"\n  confound check: r_cross vs r_replicate (reproducibility) "
          f"rho={rho_eff.statistic:+.3f}")
    print("    a strong positive here is why raw cross-context correlation is "
          "misleading:\n    it partly measures how reproducible the compound is,"
          " not how much it transfers.")

    known = v[v["moa-fine"].notna() & (v["moa-fine"] != "unclear")]
    rank = {}
    for m in ("cdi", "r_cross", "transfer"):
        g = known.groupby("moa-fine")[m].agg(["median", "size"]).query("size>=3")
        rank[m] = g["median"]
    R = pd.DataFrame(rank).dropna()
    # r_cross and transfer are inverted relative to cdi (high = conserved)
    R["r_cross_inv"] = -R.r_cross
    R["transfer_inv"] = -R.transfer
    print(f"\nmechanism ranking agreement ({len(R)} mechanisms):")
    for a, b in (("cdi", "r_cross_inv"), ("cdi", "transfer_inv"),
                 ("r_cross_inv", "transfer_inv")):
        rho = stats.spearmanr(R[a], R[b])
        print(f"  {a:13s} vs {b:13s} rho={rho.statistic:+.3f}  p={rho.pvalue:.3f}")
    print("\ntop 6 most context-specific mechanisms by each metric:")
    for m, col in (("ours (CDI)", "cdi"), ("published (1-r_cross)", "r_cross_inv"),
                   ("disattenuated (1-transfer)", "transfer_inv")):
        print(f"  {m:26s} {list(R[col].sort_values(ascending=False).index[:6])}")
    R.to_csv(TAB / "context_metric_mechanism_ranks.csv")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.4), constrained_layout=True)
    ax[0].scatter(v.cdi, -v.r_cross, s=10, alpha=0.4, color=BLUE,
                  edgecolors="none")
    rho = stats.spearmanr(v.cdi, -v.r_cross)
    ax[0].set_xlabel("our CDI"); ax[0].set_ylabel("−(cross-context r)  [published]")
    ax[0].set_title(f"A  Ours vs published (rho={rho.statistic:+.2f})",
                    loc="left", fontweight="bold", fontsize=10)

    ax[1].scatter(v.r_replicate, v.r_cross, s=10, alpha=0.4, color=ORANGE,
                  edgecolors="none")
    ax[1].plot([0, 1], [0, 1], ls="--", color="#888888", lw=1)
    ax[1].set_xlabel("reproducibility (replicate r)")
    ax[1].set_ylabel("cross-context r")
    ax[1].set_title(f"B  Why raw correlation misleads "
                    f"(rho={rho_eff.statistic:+.2f})", loc="left",
                    fontweight="bold", fontsize=10)

    ax[2].scatter(R.cdi, R.transfer_inv, s=40, color=VIOLET, edgecolors="none")
    for m, r in R.iterrows():
        ax[2].annotate(m[:18], (r.cdi, r.transfer_inv), fontsize=6,
                       xytext=(3, 2), textcoords="offset points")
    rho = stats.spearmanr(R.cdi, R.transfer_inv)
    ax[2].set_xlabel("mechanism median CDI (ours)")
    ax[2].set_ylabel("mechanism median −transfer (disattenuated)")
    ax[2].set_title(f"C  Mechanism ranks agree? (rho={rho.statistic:+.2f})",
                    loc="left", fontweight="bold", fontsize=10)
    fig.suptitle("Does the choice of context-dependence metric change the "
                 "conclusion?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "context_metric_comparison", FIG,
                    source_data={"per_drug": res, "by_mechanism": R.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
