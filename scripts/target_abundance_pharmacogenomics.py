#!/usr/bin/env python3
"""An in-vitro analogue of target-gene pharmacogenomics.

Su et al. (Nature 653:770, 2026) ran a GWAS of GLP1-receptor-agonist response in
27,885 people and found that variation in the drugs' OWN TARGET genes — a
missense variant in GLP1R for efficacy, GIPR for tirzepatide side effects —
explains part of the inter-person variability. Their conclusion: "variation in
the drug target genes contributes to inter-person variability in response".

We cannot genotype cell lines at that resolution, but the atlas offers the
direct analogue with a different source of variation. Instead of allelic
variation across 27,885 people, we have BASELINE ABUNDANCE variation of the
target across 47 cellular contexts. So we ask, per drug:

    does a line's baseline expression of the drug's own target
    predict how strongly that line responds to the drug?

  x_L  mean baseline (DMSO) log-CPM of the drug's annotated targets in line L
  y_L  response magnitude of that drug in line L (norm of the plate-matched
       delta over responsive genes, averaged across doses)
  rho  Spearman(x, y) across lines, per drug

The null is the same statistic computed with random gene sets matched on
expression level and set size, since highly expressed genes co-vary with
library composition and would otherwise manufacture a correlation.

This also tests a second idea from that paper: they found a GLP1R x GIPR
interaction (14.8-fold odds). Our analogue is whether multi-target drugs are
better predicted by their target set than single-target drugs.

Outputs: results/tables/target_abundance_pgx.csv
         figure bundle results/figures/10_drugs/target_abundance_pgx/
"""
import argparse
import ast
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
META = ROOT / "data" / "metadata" / "metadata"
FIG = ROOT / "results" / "figures" / "10_drugs"
TAB = ROOT / "results" / "tables"
MIN_LINES = 20
N_NULL = 40
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"


def parse_targets(x):
    if isinstance(x, (list, np.ndarray)):
        return [str(t).strip().upper() for t in x if str(t).strip()]
    if not isinstance(x, str) or not x.strip() or x.strip().lower() == "none":
        return []
    try:
        v = ast.literal_eval(x)
        if isinstance(v, (list, tuple)):
            return [str(t).strip().upper() for t in v if str(t).strip()]
    except Exception:
        pass
    return [t.strip().upper() for t in x.replace(";", ",").split(",") if t.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    args = ap.parse_args()
    rng = np.random.default_rng(0)

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond)
    resp = responsive_genes(DELTA, np.ones(len(G), bool))
    mag = np.linalg.norm(DELTA[:, resp], axis=1)          # response magnitude
    del DELTA

    genes = pd.read_csv(Path(ROOT / args.pb_dir) / "genes.csv")
    dmso = cond[(cond.drug == "DMSO_TF") & (cond.plate != "plate14")]
    lines = sorted(dmso.cell_line_id.unique())
    base = np.stack([np.average(X[g.row.to_numpy()], axis=0, weights=g.n_cells)
                     for _, g in dmso.groupby("cell_line_id", observed=True)])
    cpm = np.log1p(base / (base.sum(1, keepdims=True) + 1e-9) * 1e6)
    del X
    sym_pos = {}
    for i, s in enumerate(genes.gene_symbol.astype(str).str.upper()):
        sym_pos.setdefault(s, i)
    expressed = np.where(cpm.mean(0) > 0.1)[0]
    lvl = cpm.mean(0)
    print(f"{len(lines)} lines, {len(expressed)} expressed genes", flush=True)

    dm = pd.read_parquet(META / "drug_metadata.parquet")
    tg = dict(zip(dm.drug, dm.targets.map(parse_targets)))

    # per (drug, line) response magnitude, averaged over doses
    md = (pd.DataFrame({"drug": G.drug, "line": G.cell_line_id, "m": mag})
          .groupby(["drug", "line"], observed=True).m.mean())

    recs = []
    for drug, sub in md.groupby(level=0):
        ts = [sym_pos[t] for t in tg.get(drug, []) if t in sym_pos]
        s = sub.droplevel(0)
        common = [l for l in lines if l in s.index]
        if not ts or len(common) < MIN_LINES:
            continue
        y = s.reindex(common).to_numpy()
        li = [lines.index(l) for l in common]
        x = cpm[np.ix_(li, ts)].mean(1)
        if np.std(x) < 1e-9:
            continue
        rho = stats.spearmanr(x, y).statistic
        # expression-matched null
        nulls = []
        for _ in range(N_NULL):
            pick = [expressed[np.argmin(np.abs(lvl[expressed] - lvl[t]
                                               + rng.normal(0, 0.05)))]
                    for t in ts]
            xn = cpm[np.ix_(li, pick)].mean(1)
            if np.std(xn) > 1e-9:
                nulls.append(stats.spearmanr(xn, y).statistic)
        if len(nulls) < 10:
            continue
        z = (rho - np.mean(nulls)) / (np.std(nulls) + 1e-9)
        recs.append({"drug": drug, "n_targets": len(ts), "n_lines": len(common),
                     "rho": float(rho), "null_mean": float(np.mean(nulls)),
                     "z": float(z)})
    res = pd.DataFrame(recs)
    res = res.merge(dm[["drug", "moa-fine"]], on="drug", how="left")
    res.to_csv(TAB / "target_abundance_pgx.csv", index=False)

    t = stats.wilcoxon(res.rho, res.null_mean)
    print(f"\n{len(res)} drugs tested")
    print(f"median rho(target abundance, response) = {res.rho.median():+.3f} "
          f"vs null {res.null_mean.median():+.3f}")
    print(f"Wilcoxon observed vs expression-matched null: p={t.pvalue:.2e}")
    print(f"drugs with |z|>2: {(res.z.abs() > 2).sum()} "
          f"({(res.z > 2).sum()} positive, {(res.z < -2).sum()} negative)")
    print("\nstrongest target-abundance predictors of response:")
    print(res.reindex(res.z.abs().sort_values(ascending=False).index)
          .head(12)[["drug", "moa-fine", "n_targets", "rho", "z"]]
          .to_string(index=False))
    multi = res[res.n_targets > 1]
    single = res[res.n_targets == 1]
    if len(multi) > 5 and len(single) > 5:
        u = stats.mannwhitneyu(multi.z.abs(), single.z.abs())
        print(f"\nmulti-target |z| {multi.z.abs().median():.2f} vs "
              f"single-target {single.z.abs().median():.2f} (p={u.pvalue:.3f})")

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)
    axes[0].hist(res.rho, bins=40, color=BLUE, alpha=0.8, label="observed")
    axes[0].hist(res.null_mean, bins=40, color="#bdbdbd", alpha=0.6,
                 label="expression-matched null")
    axes[0].axvline(0, color="#333333", ls="--", lw=1)
    axes[0].set_xlabel("Spearman rho (target abundance vs response magnitude)")
    axes[0].set_ylabel("drugs")
    axes[0].set_title("A  Does target abundance predict response?", loc="left",
                      fontweight="bold", fontsize=10)
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].hist(res.z, bins=40, color=AQUA, alpha=0.85)
    for v in (-2, 2):
        axes[1].axvline(v, color="#888888", ls="--", lw=0.9)
    axes[1].set_xlabel("z vs expression-matched null")
    axes[1].set_ylabel("drugs")
    axes[1].set_title("B  Per-drug significance", loc="left",
                      fontweight="bold", fontsize=10)

    top = res.reindex(res.z.abs().sort_values(ascending=False).index).head(12)
    yy = np.arange(len(top))[::-1]
    axes[2].barh(yy, top.z, color=[BLUE if v > 0 else ORANGE for v in top.z],
                 height=0.62)
    axes[2].set_yticks(yy, [f"{d[:22]}" for d in top.drug], fontsize=7.5)
    axes[2].set_xlabel("z")
    axes[2].set_title("C  Strongest target-abundance effects", loc="left",
                      fontweight="bold", fontsize=10)

    fig.suptitle("In-vitro analogue of target-gene pharmacogenomics: does "
                 "target abundance predict drug response across contexts?",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "target_abundance_pgx", FIG, source_data=res,
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
