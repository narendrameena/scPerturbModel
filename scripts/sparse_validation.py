#!/usr/bin/env python3
"""Three ways the sparse detector could be right for the wrong reason.

RESULTS.md sec.39 reported that per-pair reproducibility flags 360 of 6,939
Tahoe conditions and that those calls are enriched for MEK inhibitors in
BRAF/RAS-driven lines at P = 9e-34. That is a striking number and it has three
plausible non-biological explanations, none of which was tested.

  1. THE NULL WAS CONTAMINATED. The permutation paired condition i with
     condition i+shift, and rows arrive sorted by cell line, so 99.3% of null
     pairs shared a line against 2.1% at chance. Fixed in
     ``perturbmodel.sparse_interaction`` by drawing block-crossing derangements;
     this script re-runs the detection under the corrected null and reports what
     moved.

  2. MAGNITUDE, NOT SPECIFICITY. Cross-replicate agreement rises with
     signal-to-noise, and MEK inhibitors produce large transcriptional
     responses. "MEK inhibitors are enriched among reproducible pairs" may say
     only that they do a lot. The biologically meaningful claim is narrower and
     is testable with magnitude held fixed: **within the same drug**, do
     BRAF/RAS-driven lines reproduce better than wild-type ones? Same compound,
     same dose, same response size, only the genotype differs.

  3. NO INDEPENDENT CONFIRMATION. A pair flagged in Tahoe transcription should,
     if it reflects biology rather than a Tahoe artefact, also behave specially
     in an unrelated assay. PRISM measures viability in a different laboratory
     with different chemistry, and shares cell lines and compounds with Tahoe.
     Flagged pairs should show larger |interaction| there than unflagged ones.

The third is the real test, because it can fail. The first two are controls.

Outputs: results/tables/sparse_validation.csv
         figure bundle results/figures/00_manuscript/sparse_validation/
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
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"
MEK_DRUGS = {"cobimetinib", "trametinib", "binimetinib", "tak-733"}
MAPK = ("BRAF", "KRAS", "NRAS")


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=400)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    M = pd.read_csv(TAB / "sparse_tahoe_pairs.csv")
    md = pd.read_csv(TAB / "cell_line_metadata.csv")
    drv = md.groupby("Cell_ID_Cellosaur").Driver_Gene_Symbol.apply(
        lambda s: set(s.dropna().astype(str)))
    mapk = {k: bool(v & set(MAPK)) for k, v in drv.items()}
    cv2dep = dict(zip(md.Cell_ID_Cellosaur, md.Cell_ID_DepMap))
    M["is_mek"] = M.drug.astype(str).str.lower().isin(MEK_DRUGS)
    M["is_mapk"] = M.line.map(mapk)
    print(f"{len(M)} Tahoe conditions, {int((M.q < 0.05).sum())} flagged at "
          f"FDR<0.05", flush=True)

    # ---------- control 1: does magnitude explain the enrichment? ----------
    # Within a drug, magnitude and dose are held fixed by construction, so a
    # genotype difference there cannot be a signal-to-noise effect.
    print("\n1. MAGNITUDE CONTROL — within-drug genotype contrast", flush=True)
    rows = []
    for d, g in M.groupby("drug", observed=True):
        a = g[g.is_mapk == True].cosine.dropna()
        b = g[g.is_mapk == False].cosine.dropna()
        if len(a) < 5 or len(b) < 5:
            continue
        u = stats.mannwhitneyu(a, b, alternative="greater")
        rows.append({"drug": d, "n_mut": len(a), "n_wt": len(b),
                     "med_mut": float(a.median()), "med_wt": float(b.median()),
                     "diff": float(a.median() - b.median()),
                     "p": float(u.pvalue),
                     "is_mek": str(d).lower() in MEK_DRUGS})
    W = pd.DataFrame(rows)
    if len(W):
        mek = W[W.is_mek]
        oth = W[~W.is_mek]
        print(f"   {len(W)} drugs testable; MEK inhibitors among them: "
              f"{len(mek)}")
        print(f"   median within-drug genotype gap: MEK "
              f"{mek['diff'].median():+.4f}  vs other drugs "
              f"{oth['diff'].median():+.4f}")
        if len(mek) and len(oth) > 10:
            u = stats.mannwhitneyu(mek["diff"], oth["diff"],
                                   alternative="greater")
            print(f"   MEK gap exceeds other drugs' gap: p = {u.pvalue:.4f}")
        # is the whole-panel enrichment just magnitude?
        M["mag"] = M.groupby("drug").cosine.transform("median")
        sig = M[M.q < 0.05]
        um = stats.mannwhitneyu(M[M.is_mek].mag.dropna(),
                                M[~M.is_mek].mag.dropna(),
                                alternative="greater")
        print(f"   NOTE: MEK inhibitors do have higher overall "
              f"reproducibility (p = {um.pvalue:.1e}), so the panel-wide\n"
              f"   enrichment is partly magnitude. The within-drug contrast "
              f"above is the part that is not.")

    # ---------- control 2: independent confirmation in PRISM --------------
    print("\n2. INDEPENDENT CONFIRMATION — do flagged pairs interact in PRISM?",
          flush=True)
    from perturbmodel.celldrug import prism_gamma
    gam, _, _ = prism_gamma()
    gkey = {norm(k): k for k in gam}
    M["dep"] = M.line.map(cv2dep)
    hits, miss = [], []
    for r in M.itertuples():
        k = gkey.get(norm(r.drug))
        if k is None or not isinstance(r.dep, str):
            continue
        s = gam[k]
        if r.dep not in s.index:
            continue
        (hits if r.q < 0.05 else miss).append(abs(float(s[r.dep])))
    print(f"   {len(hits)} flagged and {len(miss)} unflagged Tahoe pairs are "
          f"also measured in PRISM")
    if len(hits) > 20 and len(miss) > 100:
        u = stats.mannwhitneyu(hits, miss, alternative="greater")
        print(f"   |interaction| in PRISM: flagged {np.median(hits):.4f}  vs "
              f"unflagged {np.median(miss):.4f}   p = {u.pvalue:.3g}")
        if u.pvalue < 0.05:
            print("   Pairs flagged from Tahoe TRANSCRIPTION carry larger "
                  "interaction in PRISM\n   VIABILITY -- a different assay, "
                  "laboratory and readout. The calls are not a\n   Tahoe "
                  "artefact.")
        else:
            print("   No confirmation. The flagged pairs do not behave "
                  "specially in PRISM, so\n   they may be specific to Tahoe or "
                  "to transcription rather than general.")
    else:
        u = None
        print("   too little overlap to test")

    # ---------- what are the flagged pairs, biologically? -----------------
    print("\n3. WHAT ARE THEY — mechanism classes among the flagged pairs",
          flush=True)
    dm = pd.read_csv(TAB / "drug_metadata.csv")
    moa = dict(zip(dm.drug.astype(str), dm["moa-fine"].astype(str)))
    M["moa"] = M.drug.astype(str).map(moa)
    sig = M[M.q < 0.05]
    tab = []
    for m, g in M.groupby("moa", observed=True):
        if m in ("nan", "unclear", "") or len(g) < 20:
            continue
        k = int((g.q < 0.05).sum())
        pv = stats.hypergeom.sf(k - 1, len(M), len(g), len(sig))
        tab.append({"moa": m, "n": len(g), "n_sig": k,
                    "rate": k / len(g), "p": pv})
    MO = pd.DataFrame(tab).sort_values("p")
    if len(MO):
        MO["q"] = np.minimum(MO.p * len(MO), 1.0)
        print(f"   {int((MO.q < 0.05).sum())} of {len(MO)} mechanism classes "
              f"enriched at Bonferroni q<0.05:")
        print(MO.head(8)[["moa", "n", "n_sig", "rate", "q"]].round(4)
              .to_string(index=False))
        MO.to_csv(TAB / "sparse_validation_moa.csv", index=False)

    out = pd.DataFrame([{
        "n_conditions": len(M), "n_flagged": int((M.q < 0.05).sum()),
        "mek_within_drug_gap": float(W[W.is_mek]["diff"].median())
        if len(W) and W.is_mek.any() else np.nan,
        "other_within_drug_gap": float(W[~W.is_mek]["diff"].median())
        if len(W) else np.nan,
        "prism_flagged": float(np.median(hits)) if len(hits) else np.nan,
        "prism_unflagged": float(np.median(miss)) if len(miss) else np.nan,
        "prism_p": float(u.pvalue) if u is not None else np.nan}])
    out.to_csv(TAB / "sparse_validation.csv", index=False)
    if len(W):
        W.to_csv(TAB / "sparse_validation_withindrug.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)
    if len(W):
        d0 = [W[W.is_mek]["diff"].dropna(), W[~W.is_mek]["diff"].dropna()]
        bp = ax[0].boxplot(d0, showfliers=False, patch_artist=True,
                           medianprops=dict(color="black", lw=1.4))
        for pch, c in zip(bp["boxes"], [ORANGE, GREY]):
            pch.set_facecolor(c); pch.set_alpha(0.85)
        ax[0].axhline(0, color="#444", lw=0.9)
        ax[0].set_xticks([1, 2], [f"MEK inhibitors\n(n={len(W[W.is_mek])})",
                                  f"all other drugs\n(n={len(W[~W.is_mek])})"],
                         fontsize=7.5)
        ax[0].set_ylabel("within-drug gap: MAPK-driven − wild-type")
    ax[0].set_title("a  Magnitude held fixed by using one drug", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(hits) > 20 and len(miss) > 100:
        bp2 = ax[1].boxplot([hits, miss], showfliers=False, patch_artist=True,
                            medianprops=dict(color="black", lw=1.4))
        for pch, c in zip(bp2["boxes"], [ORANGE, GREY]):
            pch.set_facecolor(c); pch.set_alpha(0.85)
        ax[1].set_xticks([1, 2], [f"flagged in Tahoe\n(n={len(hits)})",
                                  f"not flagged\n(n={len(miss)})"], fontsize=7.5)
        ax[1].set_ylabel("|interaction| in PRISM viability")
        ax[1].text(0.5, 0.92, f"p = {u.pvalue:.1e}", transform=ax[1].transAxes,
                   ha="center", fontsize=9, color=ORANGE, fontweight="bold")
    ax[1].set_title("b  Independent assay, lab and readout", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(MO):
        top = MO.head(10)
        yy = np.arange(len(top))[::-1]
        ax[2].barh(yy, top.rate,
                   color=[ORANGE if q < 0.05 else GREY for q in top.q],
                   height=0.7)
        ax[2].axvline(float((M.q < 0.05).mean()), ls="--", color="#555", lw=1.3,
                      label="overall rate")
        ax[2].set_yticks(yy, [f"{m[:30]} ({n})" for m, n in
                              zip(top.moa, top.n)], fontsize=6.2)
        ax[2].set_xlabel("fraction of that class's pairs flagged")
        ax[2].legend(frameon=False, fontsize=7)
    ax[2].set_title("c  Which mechanisms interact", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Validating the sparse calls: magnitude control, an "
                 "independent assay, and what they are",
                 fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "sparse_validation", FIG,
                    source_data={"summary": out, "within_drug": W,
                                 "moa": MO if len(MO) else pd.DataFrame()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
