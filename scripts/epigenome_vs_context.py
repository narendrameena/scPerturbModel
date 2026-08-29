#!/usr/bin/env python3
"""Does a cell line's EPIGENETIC state predict how context-specific its drug responses are?

Three independent analyses have now failed to explain the line x drug
interaction from genetics: driver-mutation scans (0/825 at FDR), target
mutation frequency, and target abundance. The response programs themselves are
regulatory (epithelial-vs-neuronal identity, EMT, chromatin/senescence), which
predicts that the operative variable is chromatin state rather than genotype.

ENCODE has ATAC for only 3 atlas lines and DepMap blocks scripted download, but
CCLE bisulfite profiles of 928 lines are on figshare and cover 45 of our 50.
Critically these include not just methylation LEVEL but epigenetic
HETEROGENEITY metrics — proportion of discordant reads (PDR), methylation
entropy (ME), methylation haplotype load (MHL), local pairwise methylation
discordance (LPMD) — i.e. how disordered the chromatin state is within a
population, which is exactly the quantity a "regulatory plasticity" hypothesis
should implicate.

Per line we compute a context-deviance index

    line_CDI = reproducible context variance / (additive variance + context)

using cross-plate dose pairs so plate structure and noise cancel, then test it
against each epigenetic feature across lines (Spearman, BH-corrected).
Sequencing depth per line is included as a control, since a poorly covered line
would show inflated apparent context deviance for purely technical reasons.

Outputs: results/tables/epigenome_vs_context.csv
         figure bundle results/figures/11_epigenome/epigenome_vs_context/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.evaluation.delta_eval import (additive_prior, build_deltas,
                                                load_pseudobulk,
                                                responsive_genes)
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
EPI = ROOT / "data" / "external" / "ccle_epigenome"
FIG = ROOT / "results" / "figures" / "11_epigenome"
TAB = ROOT / "results" / "tables"
METRICS = {"beta": "methylation level", "pdr": "discordant reads (heterogeneity)",
           "me": "methylation entropy", "mhl": "haplotype load",
           "lpmd": "local pairwise discordance"}
FEATURES = ["genomewide", "promoter", "cgi", "cpg_shore"]
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def norm(s):
    return str(s).upper().replace("-", "").replace(" ", "").replace(".", "").replace("/", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond, keep_plate=True)
    resp = responsive_genes(DELTA, np.ones(len(G), bool))
    P = additive_prior(G, DELTA, np.ones(len(G), bool), loo=True)[:, resp]
    R = (DELTA[:, resp] - P).astype(np.float32)
    Pv = P.astype(np.float32)
    del DELTA, P
    print(f"{len(G)} conditions, {len(resp)} genes", flush=True)

    idx = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                        "drug": G.drug, "plate": G.plate})
    recs = []
    for ln, gl in idx.groupby("line", observed=True):
        cov, npair = 0.0, 0
        for _, gd in gl.groupby("drug", observed=True):
            v = gd.i.to_numpy(); pl = gd.plate.to_numpy()
            for a in range(len(v)):
                for b in range(a + 1, len(v)):
                    if pl[a] == pl[b]:
                        continue          # exclude shared plate structure
                    cov += float(np.mean(R[v[a]] * R[v[b]])); npair += 1
        if npair < 30:
            continue
        context = max(cov / npair, 0.0)
        additive = float(np.mean(Pv[gl.i.to_numpy()] ** 2))
        recs.append({"cell_line_id": ln, "n_conditions": len(gl),
                     "additive_var": additive, "context_var": context,
                     "line_cdi": context / (additive + context)
                     if additive + context > 0 else np.nan})
    cdi = pd.DataFrame(recs).dropna(subset=["line_cdi"])

    cpc = pd.read_csv(ROOT / "results/tables/cells_per_condition.csv")
    nm = cpc.drop_duplicates("cell_line").set_index("cell_line").cell_name
    depth = cpc.groupby("cell_line").n_cells.sum()
    cdi["cell_name"] = cdi.cell_line_id.map(nm)
    cdi["total_cells"] = cdi.cell_line_id.map(depth)
    cdi["key"] = cdi.cell_name.map(norm)
    print(f"{len(cdi)} lines with a context-deviance index; "
          f"median {cdi.line_cdi.median():.3f}")

    # ---------------- epigenetic features ----------------
    feats = {}
    for m in METRICS:
        f = EPI / f"ccle.{m}.csv"
        if not f.exists():
            continue
        d = pd.read_csv(f)
        d["key"] = d.cell_line_name.map(norm)
        d = d.drop_duplicates("key").set_index("key")
        for col in FEATURES:
            if col in d.columns:
                feats[f"{m}_{col}"] = pd.to_numeric(d[col], errors="coerce")
    E = pd.DataFrame(feats)
    merged = cdi.set_index("key").join(E, how="inner").dropna(subset=["line_cdi"])
    print(f"{len(merged)} lines with both context index and epigenome\n")

    rows = []
    for c in E.columns:
        v = merged[[c, "line_cdi"]].dropna()
        if len(v) < 20:
            continue
        rho = stats.spearmanr(v[c], v.line_cdi)
        # partial control: does it survive adjusting for sequencing depth?
        vv = merged[[c, "line_cdi", "total_cells"]].dropna()
        rank = vv.rank()
        resid = {}
        for col in (c, "line_cdi"):
            b = np.polyfit(rank.total_cells, rank[col], 1)
            resid[col] = rank[col] - np.polyval(b, rank.total_cells)
        rho_adj = stats.spearmanr(resid[c], resid["line_cdi"])
        rows.append({"feature": c, "metric": c.split("_")[0],
                     "region": "_".join(c.split("_")[1:]), "n": len(v),
                     "rho": rho.statistic, "p": rho.pvalue,
                     "rho_depth_adj": rho_adj.statistic, "p_depth_adj": rho_adj.pvalue})
    res = pd.DataFrame(rows)
    if len(res):
        p = res.p.to_numpy(); n = len(p); q = np.empty(n); prev = 1.0
        for r_, i in enumerate(np.argsort(p)[::-1]):
            prev = min(prev, p[i] * n / (n - r_)); q[i] = prev
        res["q"] = q
        res = res.sort_values("p")
        res.to_csv(TAB / "epigenome_vs_context.csv", index=False)
        print("epigenetic state vs context-deviance (Spearman across lines):")
        print(res.head(14).round(4).to_string(index=False))
        print(f"\n{(res.q < 0.1).sum()} of {len(res)} features at FDR<0.10")
        dep = stats.spearmanr(merged.total_cells, merged.line_cdi)
        print(f"control - sequencing depth vs context index: "
              f"rho={dep.statistic:+.3f}, p={dep.pvalue:.3f}")

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)
    top = res.head(12).iloc[::-1]
    yy = np.arange(len(top))
    axes[0].barh(yy, top.rho, color=[BLUE if v > 0 else ORANGE for v in top.rho],
                 height=0.62)
    axes[0].set_yticks(yy, [f"{r.metric}·{r.region}" for r in top.itertuples()],
                       fontsize=7.5)
    axes[0].set_xlabel("Spearman rho vs context-deviance")
    axes[0].set_title("A  Epigenetic state vs context-dependence", loc="left",
                      fontweight="bold", fontsize=10)

    best = res.iloc[0]
    v = merged[[best.feature, "line_cdi"]].dropna()
    axes[1].scatter(v[best.feature], v.line_cdi, s=26, alpha=0.75, color=AQUA,
                    edgecolors="none")
    axes[1].set_xlabel(f"{best.metric} · {best.region}")
    axes[1].set_ylabel("line context-deviance index")
    axes[1].set_title(f"B  Strongest: rho={best.rho:+.2f}, p={best.p:.1e}",
                      loc="left", fontweight="bold", fontsize=10)

    axes[2].scatter(merged.total_cells, merged.line_cdi, s=26, alpha=0.7,
                    color="#9e9e9e", edgecolors="none")
    axes[2].set_xscale("log")
    axes[2].set_xlabel("total cells profiled for the line (control)")
    axes[2].set_ylabel("line context-deviance index")
    dep = stats.spearmanr(merged.total_cells, merged.line_cdi)
    axes[2].set_title(f"C  Not a coverage artefact (rho={dep.statistic:+.2f})",
                      loc="left", fontweight="bold", fontsize=10)

    fig.suptitle("Is context-dependent drug response written in the epigenome?",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "epigenome_vs_context", FIG,
                    source_data={"per_feature": res,
                                 "per_line": merged.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
