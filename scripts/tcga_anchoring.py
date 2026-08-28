#!/usr/bin/env python3
"""Do the cell-line context-response programs have a footprint in patient tumours?

The components were defined by how cell lines differ from one another in their
drug response. If they are real transcriptional programs rather than cell-line
culture artefacts, then in primary tumours they should (a) be coherent — their
genes should covary — and (b) vary between patients. Two further questions
follow: are they associated with driver genotype in a cohort ~200x larger than
our cell-line panel, and with survival?

Scoring: each component is reduced to its top +/- loading genes and scored per
tumour as mean(z of positive genes) - mean(z of negative genes), with z taken
WITHIN cancer type so lineage differences cannot create the signal.

Tests
  coherence   mean pairwise |r| among the component's genes across tumours,
              against a null of size-matched random gene sets
  variability SD of the within-cancer-type score
  genotype    Mann-Whitney of score by non-silent mutation status per driver,
              pooled across cancer types (BH-corrected)
  survival    log-rank on median-split scores within each cancer type,
              combined across types by Stouffer

Outputs: results/tables/tcga_anchoring_*.csv
         figure bundle results/figures/09_patients/tcga_anchoring/
"""
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
TCGA = ROOT / "data" / "external" / "tcga"
FIG = ROOT / "results" / "figures" / "09_patients"
TAB = ROOT / "results" / "tables"
TOP_N = 100
MIN_SAMPLES = 50
MIN_MUT = 15
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def logrank(t1, e1, t2, e2):
    """Two-sample log-rank statistic -> (chi2, p)."""
    t = np.concatenate([t1, t2]); e = np.concatenate([e1, e2])
    g = np.concatenate([np.zeros(len(t1)), np.ones(len(t2))])
    order = np.argsort(t)
    t, e, g = t[order], e[order], g[order]
    n1 = len(t1); n2 = len(t2)
    O1 = E1 = V = 0.0
    at1, at2 = n1, n2
    i = 0
    while i < len(t):
        j = i
        while j < len(t) and t[j] == t[i]:
            j += 1
        d = e[i:j].sum()
        d1 = e[i:j][g[i:j] == 0].sum()
        n = at1 + at2
        if d > 0 and n > 1:
            O1 += d1
            E1 += d * at1 / n
            V += d * (at1 / n) * (1 - at1 / n) * (n - d) / (n - 1)
        at1 -= (g[i:j] == 0).sum()
        at2 -= (g[i:j] == 1).sum()
        i = j
    if V <= 0:
        return 0.0, 1.0
    chi2 = (O1 - E1) ** 2 / V
    return float(chi2), float(stats.chi2.sf(chi2, 1))


def bh(p):
    p = np.asarray(p, float); n = len(p); q = np.empty(n); prev = 1.0
    for r, i in enumerate(np.argsort(p)[::-1]):
        prev = min(prev, p[i] * n / (n - r)); q[i] = prev
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="dev47")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    load = pd.read_csv(TAB / f"rank1_component_loadings_{args.tag}.csv")
    stats_df = pd.read_csv(TAB / f"rank1_viability_{args.tag}.csv")
    good = stats_df[stats_df.split_half_r > 0.2].component.tolist()
    comps = [f"comp{c}" for c in good]
    load["gene_symbol"] = load.gene_symbol.str.upper()

    print("loading TCGA expression ...", flush=True)
    expr = pd.read_csv(TCGA / "pancan_expr.tsv.gz", sep="\t", index_col=0)
    expr.index = expr.index.str.upper()
    expr = expr[~expr.index.duplicated()]
    surv = pd.read_csv(TCGA / "survival.tsv", sep="\t", index_col=0)
    common = [s for s in expr.columns if s in surv.index]
    expr = expr[common]
    surv = surv.loc[common]
    ctype = surv["cancer type abbreviation"]
    print(f"{expr.shape[0]} genes x {expr.shape[1]} tumours, "
          f"{ctype.nunique()} cancer types", flush=True)

    # z-score within cancer type so lineage cannot generate the signal
    Z = pd.DataFrame(np.nan, index=expr.index, columns=expr.columns,
                     dtype=np.float32)
    for ct, idx in ctype.groupby(ctype).groups.items():
        idx = [i for i in idx if i in expr.columns]
        if len(idx) < MIN_SAMPLES:
            continue
        sub = expr[idx]
        Z[idx] = sub.sub(sub.mean(axis=1), axis=0).div(
            sub.std(axis=1).replace(0, np.nan), axis=0).astype(np.float32)
    Z = Z.dropna(axis=1, how="all")
    keep_ct = ctype.loc[Z.columns]
    print(f"{Z.shape[1]} tumours retained after per-cancer-type z-scoring")

    # ---------------- score tumours, test coherence ----------------
    scores, coh_recs = {}, []
    universe = [g for g in load.gene_symbol if g in Z.index]
    for c, col in zip(good, comps):
        sub = load[load.gene_symbol.isin(Z.index)].sort_values(col)
        neg = sub.gene_symbol.head(TOP_N).tolist()
        pos = sub.gene_symbol.tail(TOP_N).tolist()
        scores[f"C{c}"] = (Z.loc[pos].mean() - Z.loc[neg].mean()).astype(float)
        # coherence: mean |r| among component genes vs random gene sets
        gs = pos + neg
        sample_cols = rng.choice(Z.columns, min(2000, Z.shape[1]), replace=False)
        M = Z.loc[gs, sample_cols].to_numpy()
        R = np.corrcoef(M)
        obs = float(np.nanmean(np.abs(R[np.triu_indices_from(R, 1)])))
        null = []
        for _ in range(20):
            rg = rng.choice(universe, len(gs), replace=False)
            Rn = np.corrcoef(Z.loc[rg, sample_cols].to_numpy())
            null.append(np.nanmean(np.abs(Rn[np.triu_indices_from(Rn, 1)])))
        coh_recs.append({"component": f"C{c}", "coherence": obs,
                         "null_mean": float(np.mean(null)),
                         "null_sd": float(np.std(null)),
                         "z": float((obs - np.mean(null)) / (np.std(null) + 1e-9))})
    S = pd.DataFrame(scores)
    coh = pd.DataFrame(coh_recs)
    print("\ncoherence in patient tumours (mean |r| among component genes):")
    print(coh.round(3).to_string(index=False))
    S.assign(cancer_type=keep_ct).to_csv(TAB / "tcga_component_scores.csv")
    coh.to_csv(TAB / "tcga_anchoring_coherence.csv", index=False)

    # ---------------- genotype association ----------------
    mut = pd.read_csv(TCGA / "mutation.tsv.gz", sep="\t", index_col=0)
    mut = mut[[c for c in mut.columns if c in S.index]]
    g_recs = []
    for gene in mut.index:
        row = mut.loc[gene].dropna()
        m = row[row > 0].index
        w = row[row == 0].index
        if len(m) < MIN_MUT or len(w) < MIN_MUT:
            continue
        for comp in S.columns:
            a = S[comp].reindex(m).dropna()
            b = S[comp].reindex(w).dropna()
            if len(a) < MIN_MUT or len(b) < MIN_MUT:
                continue
            _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            g_recs.append({"gene": gene, "component": comp, "n_mut": len(a),
                           "n_wt": len(b),
                           "effect": float(a.median() - b.median()),
                           "p": float(p)})
    gdf = pd.DataFrame(g_recs)
    if len(gdf):
        gdf["q"] = bh(gdf.p.to_numpy())
        gdf = gdf.sort_values("p")
        gdf.to_csv(TAB / "tcga_anchoring_genotype.csv", index=False)
        print(f"\ngenotype association: {len(gdf)} tests, "
              f"{(gdf.q < 0.05).sum()} at FDR<0.05")
        print(gdf.head(12).round(4).to_string(index=False))

    # ---------------- survival ----------------
    s_recs = []
    os_t = pd.to_numeric(surv.loc[S.index, "OS.time"], errors="coerce")
    os_e = pd.to_numeric(surv.loc[S.index, "OS"], errors="coerce")
    for comp in S.columns:
        zs, ns = [], 0
        for ct, idx in keep_ct.groupby(keep_ct).groups.items():
            idx = [i for i in idx if i in S.index]
            v = S[comp].reindex(idx).dropna()
            t = os_t.reindex(v.index); e = os_e.reindex(v.index)
            ok = t.notna() & e.notna()
            v, t, e = v[ok], t[ok], e[ok]
            if len(v) < MIN_SAMPLES:
                continue
            hi = v > v.median()
            chi2, p = logrank(t[hi].to_numpy(), e[hi].to_numpy(),
                              t[~hi].to_numpy(), e[~hi].to_numpy())
            med_hi = t[hi][e[hi] == 1].median()
            med_lo = t[~hi][e[~hi] == 1].median()
            sign = 1.0 if (med_hi or 0) < (med_lo or 0) else -1.0
            zs.append(sign * stats.norm.isf(max(p, 1e-300) / 2))
            ns += 1
        if ns >= 5:
            z = float(np.sum(zs) / np.sqrt(len(zs)))
            s_recs.append({"component": comp, "n_cancer_types": ns,
                           "stouffer_z": z,
                           "p": float(2 * stats.norm.sf(abs(z)))})
    sdf = pd.DataFrame(s_recs)
    if len(sdf):
        sdf["q"] = bh(sdf.p.to_numpy())
        sdf.to_csv(TAB / "tcga_anchoring_survival.csv", index=False)
        print("\nsurvival (Stouffer across cancer types; z>0 = high score worse):")
        print(sdf.round(4).to_string(index=False))

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3), constrained_layout=True)

    x = np.arange(len(coh))
    axes[0].bar(x - 0.2, coh.coherence, width=0.4, color=BLUE, label="observed")
    axes[0].bar(x + 0.2, coh.null_mean, width=0.4, color="#bdbdbd",
                label="random gene sets")
    axes[0].set_xticks(x, coh.component)
    axes[0].set_ylabel("mean |r| among component genes")
    axes[0].set_title("A  Are the programs coherent in tumours?", loc="left",
                      fontweight="bold", fontsize=10)
    axes[0].legend(frameon=False, fontsize=8)

    big = keep_ct.value_counts().head(12).index
    data = [S.loc[keep_ct[keep_ct == ct].index.intersection(S.index),
                  S.columns[0]].dropna() for ct in big]
    vp = axes[1].violinplot(data, positions=range(len(big)), widths=0.75,
                            showmedians=True, showextrema=False)
    for b in vp["bodies"]:
        b.set_facecolor(AQUA); b.set_alpha(0.6); b.set_edgecolor("none")
    vp["cmedians"].set_color("#333333")
    axes[1].set_xticks(range(len(big)), big, rotation=45, ha="right", fontsize=8)
    axes[1].set_ylabel(f"{S.columns[0]} score (within-type z)")
    axes[1].set_title("B  Variation between patients", loc="left",
                      fontweight="bold", fontsize=10)

    if len(sdf):
        axes[2].barh(range(len(sdf)), sdf.stouffer_z,
                     color=[ORANGE if z > 0 else BLUE for z in sdf.stouffer_z],
                     height=0.6)
        axes[2].set_yticks(range(len(sdf)), sdf.component)
        for v in (-1.96, 1.96):
            axes[2].axvline(v, color="#888888", ls="--", lw=0.9)
        axes[2].set_xlabel("Stouffer z  (>0 = high score, worse survival)")
        axes[2].set_title("C  Association with outcome", loc="left",
                          fontweight="bold", fontsize=10)
    else:
        axes[2].axis("off")

    fig.suptitle("Do cell-line drug-response programs appear in patient "
                 "tumours? (TCGA pan-cancer)", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "tcga_anchoring", FIG,
                    source_data={"coherence": coh,
                                 "genotype": gdf if len(gdf) else pd.DataFrame(),
                                 "survival": sdf if len(sdf) else pd.DataFrame()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
