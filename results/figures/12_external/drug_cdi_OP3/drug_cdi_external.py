#!/usr/bin/env python3
"""Per-drug context-dependence in an external atlas, and cross-system comparison.

Our mechanism result — that how much of a drug's effect transfers between
contexts is set by its mechanism, with nuclear-receptor agonists least
transferable — was derived entirely in cancer cell lines. This recomputes the
context-dependence index (CDI) per compound in an external atlas using its own
replicate axis, then asks two questions:

  1. Does the mechanism ranking hold in the new system?
  2. For compounds measured in BOTH systems, do their CDI values agree?

Question 2 is the sharper test. If context-dependence is a property of the drug
rather than of the panel, a compound that is context-specific across cancer
cell lines should also be context-specific across immune cell types — despite
different cells, platform, lab and replicate design.

  conserved(d)  mean square of the compound's context-averaged response per dose
  context(d)    reproducible interaction: covariance of residuals between
                INDEPENDENT replicates of the same (context, compound)
  CDI           context / (conserved + context)

MOA labels are taken from the Tahoe drug metadata by compound-name match, so the
two systems are annotated identically.

Outputs: results/tables/drug_cdi_<tag>.csv
         figure bundle results/figures/12_external/drug_cdi_<tag>/
"""
import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata" / "metadata"
FIG = ROOT / "results" / "figures" / "12_external"
TAB = ROOT / "results" / "tables"
N_RESP = 2000
MIN_CONTEXTS = 3
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--context", required=True)
    ap.add_argument("--perturbation", required=True)
    ap.add_argument("--replicate", required=True)
    ap.add_argument("--dose", default=None)
    ap.add_argument("--control-value", required=True)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    import scipy.sparse as sp
    A = ad.read_h5ad(args.h5ad)
    o = A.obs.copy()
    X = A.X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    if X.min() >= 0 and X.max() > 50:
        X = np.log1p(X / (X.sum(1, keepdims=True) + 1e-9) * 1e6)
    ctx, pert, rep = args.context, args.perturbation, args.replicate
    for c in (ctx, pert, rep):
        o[c] = o[c].astype(str)
    dose = args.dose or "_dose"
    o[dose] = o[args.dose].astype(str) if args.dose else "1"

    ctrl = {k: g.index.to_numpy() for k, g in
            o[o[pert] == args.control_value].groupby([ctx, rep], observed=True)}
    pos = {n: i for i, n in enumerate(o.index)}
    rows, keys = [], []
    for name, r in o[o[pert] != args.control_value].iterrows():
        c = ctrl.get((r[ctx], r[rep]))
        if c is None or len(c) == 0:
            continue
        rows.append(X[pos[name]] - X[[pos[x] for x in c]].mean(0))
        keys.append((r[ctx], r[pert], r[dose], r[rep]))
    D = np.stack(rows).astype(np.float32)
    K = pd.DataFrame(keys, columns=["ctx", "pert", "dose", "rep"])
    resp = np.sort(np.argsort(D.var(0))[-N_RESP:])
    D = D[:, resp]
    print(f"{args.tag}: {len(D)} treated profiles, {K.ctx.nunique()} contexts, "
          f"{K.pert.nunique()} perturbations, {K.rep.nunique()} replicates",
          flush=True)

    recs = []
    for drug, gd in K.groupby("pert", observed=True):
        if gd.ctx.nunique() < MIN_CONTEXTS:
            continue
        conserved, resid = 0.0, {}
        for _, gc in gd.groupby("dose", observed=True):
            ii = gc.index.to_numpy()
            mu = D[ii].mean(0)
            conserved += float(np.mean(mu ** 2)) * len(ii)
            for i in ii:
                resid[i] = D[i] - mu
        conserved /= len(gd)
        cov, npair = 0.0, 0
        for _, gc in gd.groupby("ctx", observed=True):
            v = gc.index.to_numpy(); rp = K.rep[v].to_numpy()
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if rp[a] == rp[b]:
                        continue
                    cov += float(np.mean(resid[v[a]] * resid[v[b]])); npair += 1
        if npair < 5:
            continue
        context = max(cov / npair, 0.0)
        den = conserved + context
        recs.append({"drug": drug, "n_contexts": gd.ctx.nunique(),
                     "n_pairs": npair, "conserved": conserved,
                     "context": context,
                     "cdi": context / den if den > 0 else np.nan,
                     "total_effect": den})
    res = pd.DataFrame(recs).dropna(subset=["cdi"])

    dm = pd.read_parquet(META / "drug_metadata.parquet")
    dm["key"] = dm.drug.map(norm)
    res["key"] = res.drug.map(norm)
    res = res.merge(dm[["key", "moa-fine"]], on="key", how="left")
    res.to_csv(TAB / f"drug_cdi_{args.tag}.csv", index=False)
    print(f"{len(res)} compounds scored; median CDI {res.cdi.median():.3f}; "
          f"{res['moa-fine'].notna().sum()} with a Tahoe MOA label")

    known = res[res["moa-fine"].notna() & (res["moa-fine"] != "unclear")]
    grp = (known.groupby("moa-fine").cdi.agg(["median", "size"])
           .query("size >= 2").sort_values("median", ascending=False))
    if len(grp):
        print("\ncontext-dependence by mechanism (n>=2):")
        print(grp.round(3).to_string())

    # cross-system comparison on shared compounds
    tah = pd.read_csv(TAB / "drug_context_dependence.csv")
    tah["key"] = tah.drug.map(norm)
    both = res.merge(tah[["key", "cdi", "moa-fine"]], on="key",
                     suffixes=("_ext", "_tahoe")).dropna(subset=["cdi_ext",
                                                                 "cdi_tahoe"])
    print(f"\n{len(both)} compounds measured in BOTH systems")
    if len(both) >= 8:
        rho = stats.spearmanr(both.cdi_ext, both.cdi_tahoe)
        pr = stats.pearsonr(both.cdi_ext, both.cdi_tahoe)
        print(f"CDI agreement across systems: Spearman rho={rho.statistic:+.3f} "
              f"(p={rho.pvalue:.3f}), Pearson r={pr.statistic:+.3f}")
        both.to_csv(TAB / f"drug_cdi_shared_{args.tag}.csv", index=False)
        print(both.nlargest(8, "cdi_ext")[["drug", "moa-fine_tahoe",
                                           "cdi_ext", "cdi_tahoe"]]
              .to_string(index=False))

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    n_panels = 2 + (len(both) >= 8)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 4.3),
                             constrained_layout=True, squeeze=False)
    ax = axes[0]
    if len(grp):
        top = grp.head(12)
        yy = np.arange(len(top))[::-1]
        ax[0].barh(yy, top["median"], color=BLUE, height=0.65)
        ax[0].set_yticks(yy, [f"{m[:30]} (n={int(n)})"
                              for m, n in zip(top.index, top["size"])],
                         fontsize=7.5)
        ax[0].axvline(res.cdi.median(), color="#888888", ls="--", lw=0.9)
        ax[0].set_xlabel("context-dependence index (median)")
    ax[0].set_title(f"A  Mechanism ranking in {args.tag}", loc="left",
                    fontweight="bold", fontsize=10)

    ax[1].hist(res.cdi, bins=30, color=AQUA, alpha=0.85)
    ax[1].axvline(res.cdi.median(), color="#333333", ls="--", lw=1.2,
                  label=f"median {res.cdi.median():.2f}")
    ax[1].set_xlabel("context-dependence index")
    ax[1].set_ylabel("compounds")
    ax[1].set_title(f"B  CDI distribution ({len(res)} compounds)", loc="left",
                    fontweight="bold", fontsize=10)
    ax[1].legend(frameon=False, fontsize=8)

    if len(both) >= 8:
        ax[2].scatter(both.cdi_tahoe, both.cdi_ext, s=26, alpha=0.7,
                      color=ORANGE, edgecolors="none")
        rho = stats.spearmanr(both.cdi_ext, both.cdi_tahoe)
        ax[2].set_xlabel("CDI in Tahoe-100M (cancer lines)")
        ax[2].set_ylabel(f"CDI in {args.tag}")
        ax[2].set_title(f"C  Same drugs, both systems (rho={rho.statistic:+.2f},"
                        f" p={rho.pvalue:.2f})", loc="left", fontweight="bold",
                        fontsize=9.5)

    fig.suptitle(f"Is context-dependence a property of the drug? {args.tag}",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, f"drug_cdi_{args.tag}", FIG,
                    source_data={"per_drug": res,
                                 "shared": both if len(both) else pd.DataFrame()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
