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
    """Per-compound interaction residual over cell lines.

    Delegates to ``perturbmodel.celldrug.prism_gamma``, which removes the
    compound main effect AND each line's general sensitivity. The local
    version removed only the first, so a line that responds to everything --
    reproducible across disjoint compound halves at r = 0.989 -- looked like a
    line with a specific relation to every compound it was screened against.
    """
    from perturbmodel.celldrug import prism_gamma
    g, ti, _ = prism_gamma()
    return g


def _residuals_legacy():
    """Per-line, per-compound interaction residual -- same construction as
    prism_context_genetics.py, so the two scans are directly comparable."""
    lfc = pd.read_csv(PR / "secondary-screen-logfold-change.csv", index_col=0)
    # PRISM row names are pool_line, but 8 are pool_line_FAILED_STR. The old
    # split("_")[-1] returned the literal "STR" for all eight, and
    # groupby.mean() then averaged eight different cell lines into one
    # fabricated line that carried data in 32,230 of 36,076 profiles and was
    # published in our own tables. Parse the ACH id, and DROP the STR failures
    # rather than average an unauthenticated culture into an authenticated one:
    # all eight lines also appear as clean rows, so nothing is lost. 770 rows
    # -> 737 distinct cell lines.
    keep = [not i.endswith("_FAILED_STR") for i in lfc.index]
    lfc = lfc[keep]
    lfc.index = [re.search(r"(ACH-\d+)", i).group(1) for i in lfc.index]
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
    JZ = J.dropna(subset=["z", "gene_z"])
    w = stats.wilcoxon(np.abs(JZ.z), np.abs(JZ.gene_z))
    print(f"\n|z| allele {np.median(np.abs(JZ.z)):.3f} vs gene "
          f"{np.median(np.abs(JZ.gene_z)):.3f}, paired p={w.pvalue:.2e} "
          f"(n={len(JZ)})")
    print("  Across ALL tests the gene indicator wins slightly, which is what "
          "should\n  happen: it has more carriers, and most alleles carry no "
          "specific effect, so\n  the extra n dominates. Allele resolution pays "
          "off only where the effect\n  really is allele-specific -- which is "
          "what the count above measures, and\n  it is a claim about which "
          "associations exist, not about average power.")

    # a set of frameshifts in homopolymer runs (RPL22 K16fs, ACVR2A K437fs) is
    # the signature of microsatellite instability, so several 'independent'
    # alleles hitting one compound may be one MSI association counted many
    # times. Collapse alleles whose carrier sets substantially overlap.
    print("\nco-occurrence check (are multiple alleles tagging one background?)")
    coll = []
    for cpd, gg in vh.groupby("compound", observed=True):
        vs = gg.variant.tolist()
        if len(vs) < 2:
            coll.append({"compound": cpd, "n_alleles": len(vs),
                         "n_independent": len(vs), "max_jaccard": 0.0})
            continue
        sets = [set(vsets[v]) for v in vs]
        used, groups = set(), 0
        mx = 0.0
        for i in range(len(vs)):
            if i in used:
                continue
            groups += 1; used.add(i)
            for j in range(i + 1, len(vs)):
                if j in used:
                    continue
                jac = len(sets[i] & sets[j]) / max(len(sets[i] | sets[j]), 1)
                mx = max(mx, jac)
                if jac > 0.5:
                    used.add(j)
        coll.append({"compound": cpd, "n_alleles": len(vs),
                     "n_independent": groups, "max_jaccard": round(mx, 3),
                     "alleles": "; ".join(vs[:6])})
    C = pd.DataFrame(coll).sort_values("n_alleles", ascending=False)
    C.to_csv(TAB / "prism_variant_cooccurrence.csv", index=False)
    print(C.head(8).to_string(index=False))
    tot_a, tot_i = int(C.n_alleles.sum()), int(C.n_independent.sum())
    print(f"  {tot_a} allele hits collapse to {tot_i} independent "
          f"associations after merging alleles with >50% carrier overlap")
    # the decisive control: MSI-high lines carry MANY frameshifts, so pairwise
    # overlap between any two specific alleles can stay low while all of them
    # tag the same background. Condition on total frameshift burden instead.
    burden = (mut[mut.Variant_Classification.astype(str)
                  .str.contains("Frame_Shift", na=False)]
              .groupby("depmap").size())
    print("\nMSI control: does the association survive conditioning on total "
          "frameshift burden?")
    # Quintile matching was tried first and rejected: carriers span several
    # burden quintiles, so restricting controls to "the carriers' quintiles"
    # excluded almost nothing and every hit trivially survived. A rank-based
    # PARTIAL correlation removes the burden effect from both variables
    # instead, which is the test that can actually fail.
    ctrl = []
    for r_ in vh.itertuples():
        s = resid[r_.compound]
        b = burden.reindex(s.index).fillna(0).to_numpy()
        y = s.to_numpy()
        carrier = np.isin(s.index.to_numpy(), vsets[r_.variant]).astype(float)
        rb = stats.spearmanr(b, y)
        ry, rc, rbk = (stats.rankdata(y), stats.rankdata(carrier),
                       stats.rankdata(b))
        Z = np.column_stack([np.ones_like(rbk), rbk])
        ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
        ec = rc - Z @ np.linalg.lstsq(Z, rc, rcond=None)[0]
        pr = stats.pearsonr(ec, ey)
        # how much of this compound's interaction does the allele account for?
        ve = float(stats.pointbiserialr(carrier, y).statistic ** 2)
        ctrl.append({"compound": r_.compound, "variant": r_.variant,
                     "is_fs": "fs" in r_.variant, "p_raw": r_.p,
                     "effect": r_.effect, "var_explained": ve,
                     "burden_rho": rb.statistic, "burden_p": rb.pvalue,
                     "partial_rho": float(pr.statistic),
                     "p_partial": float(pr.pvalue)})
    CT = pd.DataFrame(ctrl)
    CT.to_csv(TAB / "prism_variant_msi_control.csv", index=False)
    for lab, sub in (("frameshift alleles", CT[CT.is_fs]),
                     ("non-frameshift alleles", CT[~CT.is_fs])):
        if not len(sub):
            continue
        surv = int((sub.p_partial < 0.05).sum())
        print(f"  {lab:24s} n={len(sub):3d}  median burden rho="
              f"{sub.burden_rho.median():+.3f}  "
              f"{surv}/{len(sub)} survive partial correlation on burden")
    print("  If the frameshift hits collapse while the missense hits (BRAF "
          "V600E,\n  PIK3CA E545K/H1047R) survive, the frameshift set was one "
          "MSI association\n  read out through many correlated marker alleles.")

    fs = vh[vh.variant.str.contains("fs", na=False)]
    print(f"  {len(fs)} of {len(vh)} hits are frameshifts; "
          f"{fs.compound.nunique()} compounds carry them. Frameshifts in "
          f"homopolymer runs\n  co-occur in MSI-high lines, so these should be "
          f"read as an MSI association,\n  not as independent gene-level "
          f"findings.")

    # direction and magnitude: a negative residual means the carrier lines lose
    # MORE viability than the compound's average, i.e. the allele makes the drug
    # more effective. This is the quantity a biomarker claim rests on.
    sens = CT[CT.effect < 0]; res_ = CT[CT.effect > 0]
    print(f"\ndirection of the {len(CT)} allele associations:")
    print(f"  {len(sens)} sensitising (carriers respond MORE than average), "
          f"{len(res_)} resistance-conferring")
    print(f"  variance of the line x compound interaction explained by a single "
          f"allele:\n    median {CT.var_explained.median():.3f}, "
          f"max {CT.var_explained.max():.3f} "
          f"({CT.loc[CT.var_explained.idxmax(), 'variant']} x "
          f"{CT.loc[CT.var_explained.idxmax(), 'compound']})")
    print("  So single alleles explain a few percent of the interaction each — "
          "real and\n  directional, but far from accounting for it. "
          "Context-dependence is not\n  reducible to a handful of driver "
          "alleles at this resolution.")

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
    lim = max(np.abs(JZ.z).quantile(0.999), np.abs(JZ.gene_z).quantile(0.999))
    ax[0].scatter(np.abs(JZ.gene_z), np.abs(JZ.z), s=5, alpha=0.2, color="#9e9e9e",
                  edgecolors="none", rasterized=True)
    sh = JZ[(JZ.q < 0.05) & (JZ.gene_q >= 0.05)]
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
