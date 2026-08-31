#!/usr/bin/env python3
"""Does drug context-dependence resolve to the variant, not just the gene?

The gene-level scan (§11) recovers the clinical biomarker set, but "BRAF is
mutated" is a coarse predictor: BRAF V600E activates the kinase and sensitises
to vemurafenib, while many non-V600 BRAF alleles do not, and some are kinase
dead. Collapsing them into one indicator averages opposite effects and costs
power. The same holds for TP53, where missense alleles in the DNA-binding
domain behave differently from truncations.

PRISM plus the Cellosaurus mutation table has protein-change resolution
(Protein_Change, e.g. p.V600E) across ~740 lines, so the question is testable:

  1. do individual alleles beat their own gene-level indicator?
  2. are there associations visible only at allele resolution -- a variant that
     is significant while its gene is not? Those are the ones a gene-level
     analysis, which is what every study of this kind runs, would miss.
  3. within one gene, do different alleles disagree in direction?

Question 2 is where new biology would be, so it is separated from the
rediscovery in (1). We report allele hits whose own gene-level test fails, and
flag which are already-known pharmacogenomics (V600E) versus not.

Outputs: results/tables/prism_variant_scan.csv
         results/tables/prism_variant_vs_gene.csv
         figure bundle results/figures/13_prism/prism_variant_level/
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
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "13_prism"
TAB = ROOT / "results" / "tables"
MIN_LINES = 100
MIN_VAR = 8
MIN_GENE = 15
N_CPD = 300
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def bh(p):
    p = np.asarray(p, float); n = len(p); q = np.empty(n); prev = 1.0
    for r, i in enumerate(np.argsort(p)[::-1]):
        prev = min(prev, p[i] * n / (n - r)); q[i] = prev
    return q


def mwu_z(ranks, Mmat, n_mut, n_tot):
    n_wt = n_tot - n_mut
    Rsum = Mmat @ ranks
    mu = n_mut * (n_tot + 1) / 2.0
    sd = np.sqrt(np.maximum(n_mut * n_wt * (n_tot + 1) / 12.0, 1e-9))
    return (Rsum - mu) / sd


def residuals():
    """Per-line, per-compound interaction residual -- same construction as
    prism_context_genetics.py, so the two scans are directly comparable."""
    lfc = pd.read_csv(PR / "secondary-screen-logfold-change.csv", index_col=0)
    lfc.index = [i.split("_")[-1] for i in lfc.index]
    lfc = lfc.groupby(level=0).mean()
    ti = pd.read_csv(PR / "secondary-screen-replicate-treatment-info.csv",
                     low_memory=False)
    ti = ti[ti.column_name.isin(lfc.columns) & ti.name.notna()].copy()
    ti["rep"] = ti.detection_plate.astype(str).str.extract(r"_(X\d)")[0]
    ti = ti[ti.rep.notna()]
    ti["dose_s"] = ti.dose.round(4).astype(str)
    L = lfc.to_numpy(dtype=np.float32)
    lines = np.array(lfc.index)
    colpos = {c: i for i, c in enumerate(lfc.columns)}
    keys, mats = [], []
    for (cpd, dose, rep), g in ti.groupby(["name", "dose_s", "rep"],
                                          observed=True):
        idx = [colpos[c] for c in g.column_name]
        with np.errstate(invalid="ignore"):
            v = np.nanmean(L[:, idx], axis=1) if len(idx) > 1 else L[:, idx[0]]
        keys.append((cpd, dose, rep)); mats.append(v)
    K = pd.DataFrame(keys, columns=["compound", "dose", "rep"])
    R = np.stack(mats, axis=1)
    out = {}
    for cpd, gc in K.groupby("compound", observed=True):
        rb = []
        for dose, gd in gc.groupby("dose", observed=True):
            cols = gd.index.to_numpy()
            if len(cols) < 2:
                continue
            sub = R[:, cols]
            with np.errstate(invalid="ignore"):
                lm = np.nanmean(sub, axis=1)
            valid = np.isfinite(lm)
            n = int(valid.sum())
            if n < MIN_LINES:
                continue
            loo = (lm[valid].sum() - lm[valid]) / (n - 1)
            res = np.full(sub.shape, np.nan, dtype=np.float32)
            res[valid] = sub[valid] - loo[:, None]
            with np.errstate(invalid="ignore"):
                rb.append(pd.Series(np.nanmean(res, axis=1), index=lines))
        if len(rb) >= 2:
            s = pd.concat(rb, axis=1).mean(axis=1).dropna()
            if len(s) >= MIN_LINES:
                out[cpd] = s
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-compounds", type=int, default=N_CPD)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    resid = residuals()
    print(f"{len(resid)} compounds with interaction residuals", flush=True)

    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(r"DepMap;\s*(ACH-\d+)")[0]
    cv2dep = dict(zip(smp.Cellosaurus_ID, dep))
    mut = pd.read_csv(PR / "mutation_long.tsv.gz", sep="\t", low_memory=False)
    mut["depmap"] = mut.RRID.map(cv2dep)
    mut = mut[mut.depmap.notna()]
    dmg = mut[mut.Variant_Classification.astype(str).str.contains(
        "Missense|Nonsense|Frame_Shift|Splice|In_Frame", na=False)].copy()
    dmg["pc"] = dmg.Protein_Change.astype(str).str.strip()
    var = dmg[dmg.pc.str.match(r"^p\.", na=False)].copy()
    var["vid"] = var.Gene_symbol.astype(str) + " " + var.pc
    vc = var.groupby("vid").depmap.nunique()
    vids = sorted(vc[vc >= MIN_VAR].index)
    gsets = {g: np.array(sorted(set(d.depmap)))
             for g, d in dmg.groupby("Gene_symbol")}
    vsets = {v: np.array(sorted(set(d.depmap)))
             for v, d in var[var.vid.isin(vids)].groupby("vid")}
    print(f"{len(vids)} distinct protein-level alleles seen in >={MIN_VAR} "
          f"lines, across {len({v.split()[0] for v in vids})} genes", flush=True)

    cpds = sorted(resid, key=lambda c: -len(resid[c]))[:args.n_compounds]
    rows = []
    for cpd in cpds:
        r = resid[cpd]
        idx = r.index.to_numpy(); vals = r.to_numpy()
        ranks = stats.rankdata(vals)
        n = len(idx)
        M = np.stack([np.isin(idx, vsets[v]) for v in vids])
        nm = M.sum(1)
        keep = (nm >= MIN_VAR) & ((n - nm) >= MIN_VAR)
        if not keep.any():
            continue
        wk = np.where(keep)[0]
        z = mwu_z(ranks, M[keep].astype(np.float64), nm[keep], n)
        pv = 2 * stats.norm.sf(np.abs(z))
        for j, i in enumerate(wk):
            v = vids[i]
            g = v.split()[0]
            rows.append({"compound": cpd, "variant": v, "gene": g,
                         "n_mut": int(nm[i]), "n_lines": n,
                         "effect": float(np.median(vals[M[i]])
                                         - np.median(vals[~M[i]])),
                         "z": float(z[j]), "p": float(pv[j])})
    V = pd.DataFrame(rows)
    V["q"] = bh(V.p.to_numpy())
    V = V.sort_values("p")
    V.to_csv(TAB / "prism_variant_scan.csv.gz", index=False, compression="gzip")
    vh = V[V.q < 0.05]
    print(f"\n{len(V)} (compound, allele) tests; {len(vh)} at FDR<0.05 "
          f"across {vh.compound.nunique()} compounds and "
          f"{vh.variant.nunique()} alleles")
    print(V.head(18)[["compound", "variant", "n_mut", "n_lines", "effect",
                      "p", "q"]].to_string(index=False))

    # gene-level test for exactly the same (compound, gene) pairs, so the
    # comparison is like-for-like rather than against the earlier scan
    need = V[["compound", "gene"]].drop_duplicates()
    grows = []
    for cpd, gg in need.groupby("compound", observed=True):
        r = resid[cpd]
        idx = r.index.to_numpy(); vals = r.to_numpy(); n = len(idx)
        ranks = stats.rankdata(vals)
        gl = [g for g in gg.gene if g in gsets]
        if not gl:
            continue
        M = np.stack([np.isin(idx, gsets[g]) for g in gl])
        nm = M.sum(1)
        keep = (nm >= MIN_GENE) & ((n - nm) >= MIN_GENE)
        if not keep.any():
            continue
        wk = np.where(keep)[0]
        z = mwu_z(ranks, M[keep].astype(np.float64), nm[keep], n)
        pv = 2 * stats.norm.sf(np.abs(z))
        for j, i in enumerate(wk):
            grows.append({"compound": cpd, "gene": gl[i],
                          "gene_n_mut": int(nm[i]), "gene_z": float(z[j]),
                          "gene_p": float(pv[j])})
    G = pd.DataFrame(grows)
    G["gene_q"] = bh(G.gene_p.to_numpy())
    J = V.merge(G, on=["compound", "gene"], how="inner")
    J.to_csv(TAB / "prism_variant_vs_gene.csv", index=False)

    sharper = J[(J.q < 0.05) & (J.gene_q >= 0.05)]
    both = J[(J.q < 0.05) & (J.gene_q < 0.05)]
    print(f"\nof {len(J)} allele tests with a matched gene-level test:")
    print(f"  {len(both)} significant at BOTH allele and gene level "
          f"(rediscovery, incl. the known biomarkers)")
    print(f"  **{len(sharper)} significant at allele level while the gene-level "
          f"test fails** — invisible to a gene-level scan")
    if len(sharper):
        print(sharper.head(15)[["compound", "variant", "n_mut", "gene_n_mut",
                                "effect", "q", "gene_q"]].to_string(index=False))

    # does resolving to the allele buy power in general?
    m = J.dropna(subset=["z", "gene_z"])
    w = stats.wilcoxon(np.abs(m.z), np.abs(m.gene_z))
    print(f"\n|z| allele {np.median(np.abs(m.z)):.3f} vs gene "
          f"{np.median(np.abs(m.gene_z)):.3f}, paired p={w.pvalue:.2e} "
          f"(n={len(m)})")
    print("  a higher allele |z| despite far fewer mutant lines means the gene "
          "indicator\n  is diluting real allele-specific effects, not that the "
          "allele test is noisier.")

    # within-gene disagreement: alleles of one gene pulling opposite ways
    dis = []
    for (cpd, g), gg in J.groupby(["compound", "gene"], observed=True):
        s = gg[gg.p < 0.01]
        if s.variant.nunique() >= 2 and s.effect.min() < 0 < s.effect.max():
            dis.append({"compound": cpd, "gene": g,
                        "n_alleles": int(s.variant.nunique()),
                        "min_effect": float(s.effect.min()),
                        "max_effect": float(s.effect.max())})
    D = pd.DataFrame(dis)
    print(f"\n{len(D)} (compound, gene) pairs where alleles of the same gene "
          f"act in OPPOSITE directions (p<0.01 each)")
    if len(D):
        print(D.head(10).round(3).to_string(index=False))

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    lim = max(np.abs(m.z).quantile(0.999), np.abs(m.gene_z).quantile(0.999))
    ax[0].scatter(np.abs(m.gene_z), np.abs(m.z), s=5, alpha=0.2, color="#9e9e9e",
                  edgecolors="none", rasterized=True)
    sh = m[(m.q < 0.05) & (m.gene_q >= 0.05)]
    ax[0].scatter(np.abs(sh.gene_z), np.abs(sh.z), s=18, color=ORANGE,
                  edgecolors="none", label=f"allele only (n={len(sh)})")
    ax[0].plot([0, lim], [0, lim], ls="--", color="#555555", lw=1)
    ax[0].set_xlim(0, lim); ax[0].set_ylim(0, lim)
    ax[0].set_xlabel("|z|, gene-level indicator")
    ax[0].set_ylabel("|z|, single allele")
    ax[0].legend(frameon=False, fontsize=8)
    ax[0].set_title("A  Resolving to the allele sharpens the test", loc="left",
                    fontweight="bold", fontsize=10)

    ax[1].scatter(V.effect, -np.log10(V.p), s=5, alpha=0.25, color="#9e9e9e",
                  edgecolors="none", rasterized=True)
    ax[1].scatter(vh.effect, -np.log10(vh.p), s=16, color=VIOLET,
                  edgecolors="none", label=f"FDR<0.05 (n={len(vh)})")
    for r_ in vh.head(8).itertuples():
        ax[1].annotate(f"{r_.variant}·{r_.compound[:9]}",
                       (r_.effect, -np.log10(r_.p)), fontsize=5.5,
                       xytext=(4, 2), textcoords="offset points")
    ax[1].set_xlabel("median residual, allele carriers − others")
    ax[1].set_ylabel("−log10 p")
    ax[1].legend(frameon=False, fontsize=8)
    ax[1].set_title("B  Allele × compound associations", loc="left",
                    fontweight="bold", fontsize=10)

    cnt = [len(both), len(sharper)]
    ax[2].bar([0, 1], cnt, color=[BLUE, ORANGE], width=0.55)
    for x, v in enumerate(cnt):
        ax[2].text(x, v + max(cnt) * 0.02, str(v), ha="center", fontsize=11)
    ax[2].set_xticks([0, 1], ["also significant\nat gene level\n(rediscovery)",
                              "allele only\n(gene-level scan\nwould miss)"],
                     fontsize=8)
    ax[2].set_ylabel("associations at FDR<0.05")
    ax[2].set_title("C  What allele resolution adds", loc="left",
                    fontweight="bold", fontsize=10)
    fig.suptitle("Does drug context-dependence resolve to the variant rather "
                 "than the gene?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "prism_variant_level", FIG,
                    source_data={"allele_scan": V[V.q < 0.25],
                                 "allele_vs_gene": J[J.p < 0.01],
                                 "opposite_direction": D}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
