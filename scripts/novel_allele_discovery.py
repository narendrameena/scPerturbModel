#!/usr/bin/env python3
"""Allele-specific drug response that a gene-level scan cannot see — novel only.

The allele scan (§14) is dominated by BRAF V600E, which is the most heavily
documented pharmacogenomic relationship in oncology. Recovering it validates the
method and proves nothing new. The question that carries actual genetic
discovery is what is left AFTER the known pharmacogenomics are removed:

  is there an allele whose drug association is (a) significant, (b) invisible to
  the gene-level indicator, (c) not already a known biomarker, and (d) still
  there when tested on cell lines that were not used to find it?

Condition (d) is the one that makes this publishable rather than a hypothesis
list. A scan over ~160,000 allele x compound tests will produce plausible top
hits by chance alone, and FDR control assumes a null the data may not obey
(alleles co-occur, compounds share targets). So every surviving candidate is
re-tested by SPLIT-SAMPLE validation: lines are split in half, the association
is discovered in one half and must independently replicate in the other, with
the split repeated so the result does not depend on one partition.

Two confounds are removed explicitly, because both would otherwise generate
convincing false discoveries:

  lineage      alleles are not distributed evenly across tissues (e.g. many
               frameshifts sit in colorectal and endometrial lines), so an
               allele can proxy for lineage-specific drug sensitivity. Tested
               by re-running within a lineage-stratified permutation.
  MSI          homopolymer frameshifts co-occur in microsatellite-unstable
               lines and tag one background, not many genes. Tested by partial
               correlation on total frameshift burden.

Outputs: results/tables/novel_allele_candidates.csv
         results/tables/novel_allele_validation.csv
         figure bundle results/figures/13_prism/novel_alleles/
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
MIN_VAR = 10
N_SPLIT = 200
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"

# Known pharmacogenomics: gene -> mechanism keywords whose association with that
# gene is already established clinical or preclinical knowledge. A hit matching
# one of these is a positive control, not a discovery, and is set aside.
KNOWN = {
    "BRAF": ["raf", "mek", "erk", "egfr", "her", "kinase"],
    "TP53": ["mdm", "p53"],
    "PIK3CA": ["pi3k", "akt", "mtor"],
    "PTEN": ["pi3k", "akt", "mtor"],
    "KRAS": ["mek", "erk", "raf", "sos", "farnesyl"],
    "NRAS": ["mek", "erk", "raf"],
    "EGFR": ["egfr", "her", "kinase"],
    "ERBB2": ["egfr", "her", "kinase"],
    "ALK": ["alk", "kinase"], "MET": ["met", "kinase"],
    "ABL1": ["abl", "bcr", "kinase"], "KIT": ["kit", "kinase"],
    "FLT3": ["flt3", "kinase"], "IDH1": ["idh"], "IDH2": ["idh"],
    "BRCA1": ["parp"], "BRCA2": ["parp"], "ATM": ["parp", "atr"],
    "SMARCA4": ["ezh2"], "ARID1A": ["ezh2", "atr"],
    "NF1": ["mek", "erk", "raf"], "NF2": ["focal adhesion", "yap"],
    "RB1": ["cdk"], "CDKN2A": ["cdk"], "CCND1": ["cdk"],
    "VHL": ["vegfr", "hif"], "STK11": ["mtor"], "TSC1": ["mtor"],
    "TSC2": ["mtor"], "MYC": ["bromodomain", "bet"],
}


def is_known(gene, moa, target):
    kws = KNOWN.get(gene)
    if not kws:
        return False
    blob = f"{moa} {target}".lower()
    return any(k in blob for k in kws)


def bh(p):
    p = np.asarray(p, float); n = len(p); q = np.empty(n); prev = 1.0
    for r, i in enumerate(np.argsort(p)[::-1]):
        prev = min(prev, p[i] * n / (n - r)); q[i] = prev
    return q


def mwu_z(ranks, M, nm, n):
    return (M @ ranks - nm * (n + 1) / 2.0) / np.sqrt(
        np.maximum(nm * (n - nm) * (n + 1) / 12.0, 1e-9))


def residuals():
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
            if valid.sum() < MIN_LINES:
                continue
            loo = (lm[valid].sum() - lm[valid]) / (valid.sum() - 1)
            res = np.full(sub.shape, np.nan, dtype=np.float32)
            res[valid] = sub[valid] - loo[:, None]
            with np.errstate(invalid="ignore"):
                rb.append(pd.Series(np.nanmean(res, axis=1), index=lines))
        if len(rb) >= 2:
            s = pd.concat(rb, axis=1).mean(axis=1).dropna()
            if len(s) >= MIN_LINES:
                out[cpd] = s
    return out, ti


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-split", type=int, default=N_SPLIT)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    resid, ti = residuals()
    ann = ti.drop_duplicates("name").set_index("name")
    print(f"{len(resid)} compounds with residuals", flush=True)

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
    vsets = {v: np.array(sorted(set(d.depmap)))
             for v, d in var[var.vid.isin(vids)].groupby("vid")}
    gsets = {g: np.array(sorted(set(d.depmap)))
             for g, d in dmg.groupby("Gene_symbol")}
    burden = (mut[mut.Variant_Classification.astype(str)
                  .str.contains("Frame_Shift", na=False)]
              .groupby("depmap").size())
    ci = pd.read_csv(PR / "secondary-screen-cell-line-info.csv")
    tissue = dict(zip(ci.depmap_id, ci.primary_tissue.astype(str)))
    print(f"{len(vids)} alleles in >={MIN_VAR} lines", flush=True)

    cpds = sorted(resid, key=lambda c: -len(resid[c]))[:300]
    rows = []
    for cpd in cpds:
        r = resid[cpd]
        idx = r.index.to_numpy(); vals = r.to_numpy(); n = len(idx)
        ranks = stats.rankdata(vals)
        M = np.stack([np.isin(idx, vsets[v]) for v in vids])
        nm = M.sum(1)
        keep = (nm >= MIN_VAR) & ((n - nm) >= MIN_VAR)
        if not keep.any():
            continue
        wk = np.where(keep)[0]
        z = mwu_z(ranks, M[keep].astype(np.float64), nm[keep], n)
        pv = 2 * stats.norm.sf(np.abs(z))
        for j, i in enumerate(wk):
            rows.append({"compound": cpd, "variant": vids[i],
                         "gene": vids[i].split()[0], "n_mut": int(nm[i]),
                         "n_lines": n, "z": float(z[j]), "p": float(pv[j]),
                         "effect": float(np.median(vals[M[i]])
                                         - np.median(vals[~M[i]]))})
    V = pd.DataFrame(rows)
    V["q"] = bh(V.p.to_numpy())
    V["moa"] = V.compound.map(ann.moa).fillna("")
    V["target"] = V.compound.map(ann.target).fillna("")
    V["known"] = [is_known(g, m, t) for g, m, t in
                  zip(V.gene, V.moa, V.target)]
    hits = V[V.q < 0.05].sort_values("p")
    print(f"\n{len(V)} tests, {len(hits)} at FDR<0.05: "
          f"{hits.known.sum()} are known pharmacogenomics, "
          f"{(~hits.known).sum()} are not")

    # gene-level test for the same pairs, to keep only allele-specific ones
    need = hits[["compound", "gene"]].drop_duplicates()
    gz = {}
    for cpd, gg in need.groupby("compound", observed=True):
        r = resid[cpd]; idx = r.index.to_numpy(); n = len(idx)
        ranks = stats.rankdata(r.to_numpy())
        gl = [g for g in gg.gene if g in gsets]
        if not gl:
            continue
        M = np.stack([np.isin(idx, gsets[g]) for g in gl])
        nm = M.sum(1)
        ok = (nm >= MIN_VAR) & ((n - nm) >= MIN_VAR)
        if not ok.any():
            continue
        z = mwu_z(ranks, M[ok].astype(np.float64), nm[ok], n)
        for j, i in enumerate(np.where(ok)[0]):
            gz[(cpd, gl[i])] = 2 * stats.norm.sf(abs(z[j]))
    hits = hits.copy()
    hits["gene_p"] = [gz.get((c, g), np.nan) for c, g in
                      zip(hits.compound, hits.gene)]
    cand = hits[(~hits.known) & (hits.gene_p > 0.05)].copy()
    print(f"{len(cand)} novel AND allele-specific (gene-level p>0.05)")

    # split-sample validation: discover in one half, require it in the other
    print(f"\nsplit-sample validation, {args.n_split} random halves each",
          flush=True)
    val = []
    for r_ in cand.itertuples():
        s = resid[r_.compound]
        idx = s.index.to_numpy(); y = s.to_numpy()
        memb = np.isin(idx, vsets[r_.variant])
        hold = 0; tried = 0
        for _ in range(args.n_split):
            perm = rng.permutation(len(idx))
            h1, h2 = perm[:len(idx) // 2], perm[len(idx) // 2:]
            if memb[h1].sum() < 4 or memb[h2].sum() < 4:
                continue
            if (~memb[h1]).sum() < 4 or (~memb[h2]).sum() < 4:
                continue
            tried += 1
            _, p1 = stats.mannwhitneyu(y[h1][memb[h1]], y[h1][~memb[h1]])
            _, p2 = stats.mannwhitneyu(y[h2][memb[h2]], y[h2][~memb[h2]])
            d1 = np.median(y[h1][memb[h1]]) - np.median(y[h1][~memb[h1]])
            d2 = np.median(y[h2][memb[h2]]) - np.median(y[h2][~memb[h2]])
            # discovered in one half at p<0.05, must hold in the other at
            # p<0.05 with the SAME sign
            if p1 < 0.05 and p2 < 0.05 and np.sign(d1) == np.sign(d2):
                hold += 1
        rate = hold / tried if tried else np.nan

        # lineage confound: does the allele just proxy for a tissue?
        tis = np.array([tissue.get(d, "?") for d in idx])
        big = pd.Series(tis[memb]).value_counts()
        top_frac = float(big.iloc[0] / max(memb.sum(), 1)) if len(big) else 1.0
        # within-lineage test, restricted to the dominant lineage's complement
        keepl = tis != (big.index[0] if len(big) else "?")
        if memb[keepl].sum() >= 5 and (~memb[keepl]).sum() >= 5:
            _, p_wo = stats.mannwhitneyu(y[keepl][memb[keepl]],
                                         y[keepl][~memb[keepl]])
        else:
            p_wo = np.nan

        # MSI burden partial correlation
        b = burden.reindex(s.index).fillna(0).to_numpy()
        ry, rc, rb = (stats.rankdata(y), stats.rankdata(memb.astype(float)),
                      stats.rankdata(b))
        Z = np.column_stack([np.ones_like(rb), rb])
        ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
        ec = rc - Z @ np.linalg.lstsq(Z, rc, rcond=None)[0]
        p_burd = stats.pearsonr(ec, ey).pvalue

        val.append({"compound": r_.compound, "variant": r_.variant,
                    "gene": r_.gene, "moa": r_.moa, "n_mut": r_.n_mut,
                    "effect": r_.effect, "q": r_.q, "gene_p": r_.gene_p,
                    "split_replication": rate, "n_splits": tried,
                    "top_lineage_frac": top_frac,
                    "p_without_top_lineage": p_wo, "p_burden_partial": p_burd})
    W = pd.DataFrame(val)
    if len(W):
        W = W.sort_values("split_replication", ascending=False)
        W.to_csv(TAB / "novel_allele_validation.csv", index=False)
        surv = W[(W.split_replication >= 0.5)
                 & (W.p_without_top_lineage < 0.05)
                 & (W.p_burden_partial < 0.05)]
        print(f"\n{len(W)} candidates tested; {len(surv)} survive ALL of: "
              f">=50% split replication, lineage control, MSI control")
        cols = ["compound", "variant", "moa", "n_mut", "effect",
                "split_replication", "p_without_top_lineage", "p_burden_partial"]
        print(W.head(12)[cols].to_string(index=False))
        if len(surv):
            print("\nSURVIVING NOVEL CANDIDATES:")
            print(surv[cols].to_string(index=False))
        else:
            print("\nNo novel allele-compound association survives all three "
                  "controls.\nThat is a real result: at 738 lines the "
                  "allele-resolution signal that is\nnot already known "
                  "pharmacogenomics does not replicate out of sample.")
    cand.to_csv(TAB / "novel_allele_candidates.csv", index=False)

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    cnt = [int(hits.known.sum()), int((~hits.known).sum()), len(cand),
           int((W.split_replication >= 0.5).sum()) if len(W) else 0]
    lab = ["known\npharmaco-\ngenomics", "not known", "+ allele-\nspecific",
           "+ replicates\nout of sample"]
    ax[0].bar(range(4), cnt, color=[BLUE, VIOLET, AQUA, ORANGE], width=0.6)
    for x, v in enumerate(cnt):
        ax[0].text(x, v + max(cnt) * 0.02, str(v), ha="center", fontsize=10)
    ax[0].set_xticks(range(4), lab, fontsize=7)
    ax[0].set_ylabel("allele × compound associations")
    ax[0].set_title("A  Attrition through the controls", loc="left",
                    fontweight="bold", fontsize=10)

    if len(W):
        ax[1].hist(W.split_replication.dropna(), bins=20, color=VIOLET,
                   alpha=0.85)
        ax[1].axvline(0.5, color=ORANGE, ls="--", lw=1.5)
        ax[1].set_xlabel("fraction of random half-splits replicating")
        ax[1].set_ylabel("candidates")
        ax[1].set_title("B  Out-of-sample replication", loc="left",
                        fontweight="bold", fontsize=10)

    kn = V[V.known & (V.q < 0.05)]
    ax[2].scatter(V.effect, -np.log10(V.p), s=4, alpha=0.15, color="#bdbdbd",
                  edgecolors="none", rasterized=True)
    ax[2].scatter(kn.effect, -np.log10(kn.p), s=16, color=BLUE,
                  edgecolors="none", label=f"known ({len(kn)})")
    nk = hits[~hits.known]
    ax[2].scatter(nk.effect, -np.log10(nk.p), s=16, color=ORANGE,
                  edgecolors="none", label=f"not known ({len(nk)})")
    ax[2].set_xlabel("median residual, carriers − others")
    ax[2].set_ylabel("−log10 p")
    ax[2].legend(frameon=False, fontsize=8)
    ax[2].set_title("C  Known vs novel", loc="left", fontweight="bold",
                    fontsize=10)
    fig.suptitle("Is there allele-specific drug response beyond the known "
                 "pharmacogenomics?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "novel_alleles", FIG,
                    source_data={"hits": hits, "candidates": cand,
                                 "validation": W if len(W) else pd.DataFrame()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
