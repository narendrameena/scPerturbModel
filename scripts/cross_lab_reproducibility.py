#!/usr/bin/env python3
"""How much of a cell line's drug-specific response reproduces in another lab?

Ben-David et al. (2018, Nature) showed that "the same" cell line differs between
laboratories -- a median 19% of non-silent mutations appear in only one of
CCLE/GDSC, and 48 of 55 compounds active against some MCF7 strain were completely
inactive against another. They demonstrated it in one cell line, at one dose,
with viability. The question that matters for everyone using these atlases is the
scaled-up one, and it has not been answered: **of the line x compound interaction
we attribute to cell-line identity, what fraction reproduces across labs?**

That number sets a hard ceiling on every context-transfer model trained on one
atlas and applied elsewhere, including ours.

The design separates the two things that make a cross-lab correlation low:

  ceiling      within-lab, cross-experiment agreement — GDSC1 vs GDSC2 (both
               Sanger) and PRISM replicate plates (Broad). This is how well the
               measurement agrees with ITSELF, and it is not 1.
  cross-lab    PRISM (Broad, pooled-barcode viability) vs GDSC (Sanger, fitted
               IC50). Different institution, different assay, same phenotype.

Reporting cross-lab agreement alone would confuse laboratory divergence with
assay noise. The informative quantity is the **reproducible fraction**,
cross-lab r divided by the within-lab ceiling: 1.0 means labs agree as well as an
assay agrees with itself, 0 means nothing about the line survives the transfer.

The phenotype is the interaction residual used throughout this project -- the
line's deviation from the compound's average effect across lines -- not raw
potency. Raw potency correlates across labs largely because some compounds are
simply more toxic than others everywhere; removing the compound main effect
leaves exactly the context-specific part this project is about.

Outputs: results/tables/cross_lab_reproducibility.csv
         results/tables/cross_lab_summary.csv
         figure bundle results/figures/16_crosslab/cross_lab_reproducibility/
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
GD = ROOT / "data" / "external" / "gdsc"
FIG = ROOT / "results" / "figures" / "16_crosslab"
TAB = ROOT / "results" / "tables"
MIN_SHARED = 20
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def prism_residuals():
    """Per (line, compound) interaction residual, plus a replicate-split copy
    for the internal ceiling."""
    lfc = pd.read_csv(PR / "secondary-screen-logfold-change.csv", index_col=0)
    lfc.index = [i.split("_")[-1] for i in lfc.index]
    lfc = lfc.groupby(level=0).mean()
    ti = pd.read_csv(PR / "secondary-screen-replicate-treatment-info.csv",
                     low_memory=False)
    ti = ti[ti.column_name.isin(lfc.columns) & ti.name.notna()].copy()
    ti["rep"] = ti.detection_plate.astype(str).str.extract(r"_(X\d)")[0]
    ti = ti[ti.rep.notna()]
    L = lfc.to_numpy(dtype=np.float32)
    lines = np.array(lfc.index)
    colpos = {c: i for i, c in enumerate(lfc.columns)}

    def build(mask):
        out = {}
        sub = ti[mask]
        for cpd, g in sub.groupby("name", observed=True):
            idx = [colpos[c] for c in g.column_name]
            if len(idx) < 2:
                continue
            with np.errstate(invalid="ignore"):
                v = np.nanmean(L[:, idx], axis=1)
            ok = np.isfinite(v)
            if ok.sum() < MIN_SHARED:
                continue
            n = ok.sum()
            loo = (v[ok].sum() - v[ok]) / (n - 1)
            out[cpd] = pd.Series(v[ok] - loo, index=lines[ok])
        return out

    full = build(np.ones(len(ti), bool))
    half1 = build((ti.rep == "X1").to_numpy())
    half2 = build(ti.rep.isin(["X2", "X3"]).to_numpy())
    return full, half1, half2


def gdsc_residuals():
    """Z_SCORE is already the line's deviation from the drug mean. Kept per
    screening version so GDSC1 vs GDSC2 can serve as the within-lab control."""
    out = {}
    for tag, f in (("GDSC1", "GDSC1_fitted_dose_response_27Oct23.csv"),
                   ("GDSC2", "GDSC2_fitted_dose_response_27Oct23.csv")):
        g = pd.read_csv(GD / f, low_memory=False)
        g["k"] = g.DRUG_NAME.map(norm)
        out[tag] = g
    return out


def per_compound_corr(a, b, label, min_shared=MIN_SHARED):
    """Spearman per compound over shared cell lines."""
    rows = []
    for cpd in set(a) & set(b):
        x, y = a[cpd], b[cpd]
        j = x.index.intersection(y.index)
        if len(j) < min_shared:
            continue
        xv, yv = x[j].to_numpy(), y[j].to_numpy()
        if np.std(xv) == 0 or np.std(yv) == 0:
            continue
        r = stats.spearmanr(xv, yv)
        rows.append({"comparison": label, "compound": cpd, "n_lines": len(j),
                     "rho": float(r.statistic), "p": float(r.pvalue)})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-shared", type=int, default=MIN_SHARED)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    print("building PRISM residuals ...", flush=True)
    p_full, p_h1, p_h2 = prism_residuals()
    print(f"  {len(p_full)} compounds; replicate halves "
          f"{len(p_h1)}/{len(p_h2)}", flush=True)

    print("loading GDSC ...", flush=True)
    G = gdsc_residuals()
    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(r"DepMap;\s*(ACH-\d+)")[0]
    cos = {}
    for d, xr in zip(dep, smp.cross_references.astype(str)):
        if not isinstance(d, str):
            continue
        for cid in re.findall(r"Cosmic(?:-CLP)?;\s*(\d+)", xr):
            cos.setdefault(int(cid), d)

    def gdsc_map(tag):
        g = G[tag].copy()
        g["dep"] = g.COSMIC_ID.map(cos)
        g = g[g.dep.notna()]
        out = {}
        for k, gg in g.groupby("k", observed=True):
            v = gg.groupby("dep").Z_SCORE.mean()
            if len(v) >= args.min_shared:
                out[k] = v
        return out

    g1, g2 = gdsc_map("GDSC1"), gdsc_map("GDSC2")
    gboth = {}
    for k in set(g1) | set(g2):
        s = pd.concat([g1.get(k, pd.Series(dtype=float)),
                       g2.get(k, pd.Series(dtype=float))])
        gboth[k] = s.groupby(level=0).mean()
    print(f"  GDSC1 {len(g1)} drugs, GDSC2 {len(g2)}, union {len(gboth)}; "
          f"{len(cos)} COSMIC->DepMap mappings", flush=True)

    p_full_n = {norm(k): v for k, v in p_full.items()}
    p_h1_n = {norm(k): v for k, v in p_h1.items()}
    p_h2_n = {norm(k): v for k, v in p_h2.items()}

    rows = []
    rows += per_compound_corr(p_h1_n, p_h2_n,
                              "PRISM replicate split (within-lab ceiling)")
    rows += per_compound_corr(g1, g2, "GDSC1 vs GDSC2 (within-lab ceiling)")
    rows += per_compound_corr(p_full_n, gboth, "PRISM vs GDSC (CROSS-LAB)")
    R = pd.DataFrame(rows)
    R.to_csv(TAB / "cross_lab_reproducibility.csv", index=False)

    print("\n=== per-compound Spearman of the line-specific residual ===")
    S = (R.groupby("comparison").rho
         .agg(["median", "mean", "size",
               lambda s: float((s > 0).mean())])
         .rename(columns={"<lambda_0>": "frac_positive"}))
    print(S.round(3).to_string())

    ceil_p = R[R.comparison.str.startswith("PRISM replicate")].rho.median()
    ceil_g = R[R.comparison.str.startswith("GDSC1")].rho.median()
    cross = R[R.comparison.str.contains("CROSS-LAB")].rho.median()
    ceiling = float(np.sqrt(max(ceil_p, 1e-9) * max(ceil_g, 1e-9)))
    frac = cross / ceiling if ceiling > 0 else np.nan
    print(f"\nwithin-lab ceilings: PRISM {ceil_p:.3f}, GDSC {ceil_g:.3f} "
          f"(geometric mean {ceiling:.3f})")
    print(f"cross-lab PRISM vs GDSC: {cross:.3f}")
    print(f"REPRODUCIBLE FRACTION = {frac:.1%}")
    print("  i.e. of the line-specific response an assay can reproduce with "
          "itself,\n  this much survives moving to another laboratory and assay.")

    # does agreement depend on how strong the compound is?
    xl = R[R.comparison.str.contains("CROSS-LAB")].copy()
    strength = {norm(k): float(np.mean(v.to_numpy() ** 2))
                for k, v in p_full.items()}
    xl["interaction_strength"] = xl.compound.map(strength)
    v = xl.dropna(subset=["interaction_strength"])
    if len(v) > 10:
        rho = stats.spearmanr(v.interaction_strength, v.rho)
        print(f"\ncross-lab agreement vs interaction strength: "
              f"rho={rho.statistic:+.3f}, p={rho.pvalue:.2e} (n={len(v)})")
        q = pd.qcut(v.interaction_strength, 4, labels=["Q1 weak", "Q2", "Q3",
                                                       "Q4 strong"])
        print(v.groupby(q, observed=True).rho.agg(["median", "size"])
              .round(3).to_string())
        print("  Compounds with a stronger line-specific component reproduce "
              "better,\n  which is what a real but noisy signal looks like — "
              "and it means the\n  reproducible fraction above is a floor for "
              "the compounds worth modelling.")

    summary = pd.DataFrame([
        {"quantity": "PRISM replicate ceiling", "value": ceil_p},
        {"quantity": "GDSC1 vs GDSC2 ceiling", "value": ceil_g},
        {"quantity": "cross-lab PRISM vs GDSC", "value": cross},
        {"quantity": "reproducible fraction", "value": frac}])
    summary.to_csv(TAB / "cross_lab_summary.csv", index=False)

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    order = ["PRISM replicate split (within-lab ceiling)",
             "GDSC1 vs GDSC2 (within-lab ceiling)",
             "PRISM vs GDSC (CROSS-LAB)"]
    data = [R[R.comparison == o].rho.dropna() for o in order]
    bp = ax[0].boxplot(data, showfliers=False, patch_artist=True,
                       medianprops=dict(color="black", lw=1.6))
    for patch, c in zip(bp["boxes"], [AQUA, BLUE, ORANGE]):
        patch.set_facecolor(c); patch.set_alpha(0.75)
    ax[0].set_xticks([1, 2, 3], ["PRISM\nreplicates", "GDSC1\nvs GDSC2",
                                 "PRISM vs GDSC\n(cross-lab)"], fontsize=8)
    ax[0].axhline(0, color="#444444", lw=0.9)
    ax[0].set_ylabel("Spearman r of the line-specific residual")
    for i, d in enumerate(data):
        ax[0].text(i + 1, np.median(d) + 0.02, f"{np.median(d):.2f}",
                   ha="center", fontsize=9, fontweight="bold")
    ax[0].set_title("A  Same assay, same lab, then another lab", loc="left",
                    fontweight="bold", fontsize=10)

    for o, c in zip(order, [AQUA, BLUE, ORANGE]):
        d = R[R.comparison == o].rho.dropna()
        ax[1].hist(d, bins=30, histtype="step", lw=2, color=c,
                   label=o.split("(")[0].strip(), density=True)
    ax[1].axvline(0, color="#444444", lw=0.9)
    ax[1].set_xlabel("per-compound Spearman r")
    ax[1].set_ylabel("density"); ax[1].legend(frameon=False, fontsize=7)
    ax[1].set_title("B  Distributions", loc="left", fontweight="bold",
                    fontsize=10)

    if len(v) > 10:
        ax[2].scatter(v.interaction_strength, v.rho, s=12, alpha=0.45,
                      color=VIOLET, edgecolors="none")
        ax[2].set_xscale("log")
        ax[2].axhline(0, color="#444444", lw=0.9)
        ax[2].axhline(ceiling, color=AQUA, ls="--", lw=1.4,
                      label="within-lab ceiling")
        ax[2].set_xlabel("strength of the line-specific component (PRISM)")
        ax[2].set_ylabel("cross-lab Spearman r")
        ax[2].legend(frameon=False, fontsize=8)
        ax[2].set_title(f"C  Stronger signal transfers better\n"
                        f"rho={rho.statistic:+.2f}", loc="left",
                        fontweight="bold", fontsize=10)
    fig.suptitle("How much of a cell line's drug-specific response survives "
                 "moving to another laboratory?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "cross_lab_reproducibility", FIG,
                    source_data={"per_compound": R, "summary": summary},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
