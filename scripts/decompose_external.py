#!/usr/bin/env python3
"""Replicate the variance decomposition in an external perturbation atlas.

Our decomposition was derived entirely on Tahoe-100M, so it needs to be shown
portable. This runs the identical analysis on other datasets by mapping their
columns onto the same three roles:

    context     the thing that varies and might modulate response
                (Tahoe: cell line;  OP3: immune cell type;  sciPlex3: cell line)
    perturbation the drug
    replicate   an independent repeat of the same context x perturbation
                (Tahoe: plate;  OP3: donor;  sciPlex3: none available)

The replicate axis is what makes the estimate trustworthy: reproducible
interaction is measured as the covariance of residuals between INDEPENDENT
replicates of the same (context, perturbation), so independent noise and
batch structure both cancel. OP3's donors are biological replicates, which is
a stronger control than Tahoe's plates.

Reported, exactly as for Tahoe:
    additive     mean square of the perturbation's own context-averaged effect
    interaction  reproducible context x perturbation variance
    the three residual correlations (across replicates / across perturbations
    within a context / across contexts for a perturbation)

Outputs: results/tables/decompose_<tag>.csv  and a figure bundle.
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
FIG = ROOT / "results" / "figures" / "12_external"
TAB = ROOT / "results" / "tables"
N_RESP = 2000
MAX_PAIRS = 4000
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def corr_rows(A, B):
    A = A - A.mean(1, keepdims=True)
    B = B - B.mean(1, keepdims=True)
    d = np.linalg.norm(A, axis=1) * np.linalg.norm(B, axis=1)
    return np.where(d > 0, (A * B).sum(1) / np.where(d > 0, d, 1), np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--context", required=True)
    ap.add_argument("--perturbation", required=True)
    ap.add_argument("--replicate", default=None)
    ap.add_argument("--dose", default=None)
    ap.add_argument("--control-value", required=True,
                    help="value of --perturbation marking untreated controls")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    import anndata as ad
    import scipy.sparse as sp
    A = ad.read_h5ad(args.h5ad)
    obs = A.obs.copy()
    X = A.X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    if X.min() >= 0 and X.max() > 50:               # looks like counts
        X = np.log1p(X / (X.sum(1, keepdims=True) + 1e-9) * 1e6)
    print(f"{args.tag}: {X.shape[0]} profiles x {X.shape[1]} genes")

    ctx, pert = args.context, args.perturbation
    rep = args.replicate if args.replicate else None
    obs[ctx] = obs[ctx].astype(str)
    obs[pert] = obs[pert].astype(str)
    if rep:
        obs[rep] = obs[rep].astype(str)
    else:
        obs["_rep"] = "single"; rep = "_rep"
    dose = args.dose if args.dose else None
    if dose:
        obs[dose] = obs[dose].astype(str)
    else:
        obs["_dose"] = "1"; dose = "_dose"
    print(f"  contexts={obs[ctx].nunique()}  perturbations={obs[pert].nunique()}"
          f"  replicates={obs[rep].nunique()}  doses={obs[dose].nunique()}")

    # replicate-matched control per (context, replicate)
    is_ctrl = obs[pert] == args.control_value
    ctrl_idx = {k: g.index.to_numpy()
                for k, g in obs[is_ctrl].groupby([ctx, rep], observed=True)}
    pos = {n: i for i, n in enumerate(obs.index)}
    rows, keys = [], []
    for i, (name, r) in enumerate(obs[~is_ctrl].iterrows()):
        c = ctrl_idx.get((r[ctx], r[rep]))
        if c is None or len(c) == 0:
            continue
        rows.append(X[pos[name]] - X[[pos[x] for x in c]].mean(0))
        keys.append((r[ctx], r[pert], r[dose], r[rep]))
    if not rows:
        print("no replicate-matched controls found"); return
    D = np.stack(rows).astype(np.float32)
    K = pd.DataFrame(keys, columns=["ctx", "pert", "dose", "rep"])
    resp = np.sort(np.argsort(D.var(0))[-N_RESP:])
    D = D[:, resp]
    print(f"  {len(D)} treated profiles with matched controls")

    # additive prior: mean over contexts for each (perturbation, dose)
    P = np.zeros_like(D)
    for _, g in K.groupby(["pert", "dose"], observed=True):
        ii = g.index.to_numpy()
        for i in ii:                                  # leave-one-context-out
            other = [j for j in ii if K.ctx[j] != K.ctx[i]]
            if other:
                P[i] = D[other].mean(0)
    R = D - P

    # reproducible interaction: independent replicates of the same (ctx, pert)
    cov, npair = 0.0, 0
    for _, g in K.groupby(["ctx", "pert"], observed=True):
        v = g.index.to_numpy(); rp = K.rep[v].to_numpy()
        for a in range(len(v)):
            for b in range(a + 1, len(v)):
                if rp[a] == rp[b]:
                    continue                          # require independence
                cov += float(np.mean(R[v[a]] * R[v[b]])); npair += 1
    interaction = max(cov / max(npair, 1), 0.0)
    additive = float(np.mean(P ** 2))
    total_repro = additive + interaction
    print(f"\n  cross-replicate (ctx,pert) pairs: {npair}")
    print(f"  additive    {additive:.5f}  ({additive/total_repro:.1%} of reproducible)")
    print(f"  interaction {interaction:.5f}  ({interaction/total_repro:.1%})")

    # three residual correlations, replicate-independent pairs only
    def sample(match, differ):
        out = []
        for _, g in K.groupby(match, observed=True):
            v = g.index.to_numpy()
            if len(v) < 2:
                continue
            for _ in range(min(len(v), 60)):
                i, j = rng.choice(v, 2, replace=False)
                if K.loc[i, differ] != K.loc[j, differ] and K.rep[i] != K.rep[j]:
                    out.append((i, j))
        return np.array(out[:MAX_PAIRS]) if out else np.empty((0, 2), int)

    tests = {"across_replicates_same_ctx_pert": (["ctx", "pert"], "rep"),
             "across_perturbations_within_ctx": (["ctx"], "pert"),
             "across_contexts_same_perturbation": (["pert"], "ctx")}
    crecs = []
    for nm, (m, d) in tests.items():
        pr = sample(m, d)
        if len(pr) < 30:
            print(f"  {nm}: too few pairs"); continue
        r = corr_rows(R[pr[:, 0]], R[pr[:, 1]])
        pi = rng.integers(0, len(R), len(pr)); pj = rng.integers(0, len(R), len(pr))
        rn = corr_rows(R[pi], R[pj])
        crecs += [{"test": nm, "kind": "observed", "r": v} for v in r]
        crecs += [{"test": nm, "kind": "null", "r": v} for v in rn]
        print(f"  {nm:36s} median r = {np.nanmedian(r):+.3f} "
              f"(null {np.nanmedian(rn):+.3f}, n={len(r)})")

    summ = pd.DataFrame([{"dataset": args.tag, "n_contexts": obs[ctx].nunique(),
                          "n_perturbations": obs[pert].nunique(),
                          "n_replicates": obs[rep].nunique(),
                          "additive": additive, "interaction": interaction,
                          "interaction_share": interaction / total_repro,
                          "n_cross_replicate_pairs": npair}])
    summ.to_csv(TAB / f"decompose_{args.tag}.csv", index=False)
    cdf = pd.DataFrame(crecs)
    if len(cdf):
        cdf.to_csv(TAB / f"decompose_{args.tag}_correlations.csv", index=False)

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    sh = [additive / total_repro, interaction / total_repro]
    axes[0].bar([0, 1], sh, width=0.6, color=[BLUE, ORANGE])
    for x, v in enumerate(sh):
        axes[0].text(x, v + 0.01, f"{v:.0%}", ha="center", fontsize=11)
    axes[0].set_xticks([0, 1], ["perturbation\n(additive)",
                                "context × perturbation\n(interaction)"])
    axes[0].set_ylabel("share of reproducible response variance")
    axes[0].set_ylim(0, 1)
    axes[0].set_title(f"A  {args.tag}", loc="left", fontweight="bold", fontsize=10)

    if len(cdf):
        order = [t for t in tests if (cdf.test == t).any()]
        pos_ = 0; ticks = []; labels = []
        for t in order:
            for kind, col in (("observed", AQUA), ("null", "#bdbdbd")):
                v = cdf[(cdf.test == t) & (cdf.kind == kind)].r.dropna()
                vp = axes[1].violinplot([v], positions=[pos_], widths=0.75,
                                        showmedians=True, showextrema=False)
                vp["bodies"][0].set_facecolor(col)
                vp["bodies"][0].set_alpha(0.65 if kind == "observed" else 0.4)
                vp["bodies"][0].set_edgecolor("none")
                vp["cmedians"].set_color("#333333")
                pos_ += 1
            ticks.append(pos_ - 1.5)
            labels.append(t.replace("_", "\n")[:34])
            pos_ += 0.6
        axes[1].axhline(0, color="#888888", lw=0.9, ls="--")
        axes[1].set_xticks(ticks, labels, fontsize=7.5)
        axes[1].set_ylabel("residual correlation")
        axes[1].set_title("B  Structure of the residual", loc="left",
                          fontweight="bold", fontsize=10)
    fig.suptitle(f"Variance decomposition replicated in {args.tag}",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, f"decompose_{args.tag}", FIG,
                    source_data={"summary": summ, "correlations": cdf},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
