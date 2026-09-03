#!/usr/bin/env python3
"""Does baseline expression predict the cell-drug RELATION, or just the cell?

RESULTS.md sec.17 reports that baseline expression predicts the residual drug
response (CV R2 +0.092) while mutations do not (-0.000), and concludes that
context-dependent response is governed by transcriptional state rather than
genotype. That conclusion rests on a residual from which only the compound main
effect had been removed. A line's general sensitivity -- how it responds to
everything, driven by growth rate, seeding density and drug metabolism -- was
still inside it.

Baseline expression predicts growth rate very well. So the reported R2 may be
almost entirely the prediction of a CELL PROPERTY, dressed up as the prediction
of a cell-drug relation. Ben-David et al. (Nature 2018) show the property is
real and large: their most resistant MCF7 strains are resistant in general,
through downregulated drug-metabolism pathways.

This script splits the question in two and asks it twice, on the same lines,
the same predictors and the same folds:

    A. predict alpha, the line's general sensitivity      -- a cell property
    B. predict gamma, the interaction after alpha is out  -- a cell-drug relation

If expression predicts A well and B near zero, the claim as written is wrong
and must be restated: transcriptional state predicts how sensitive a cell is,
not which drug it will respond to. If it predicts B too, the claim survives on
a residual that can no longer be confused with general sensitivity.

The genotype blocks are carried through the identical split so that the
state-versus-genotype contrast is made where it belongs -- inside B.

Outputs: results/tables/state_vs_genotype_split.csv
         figure bundle results/figures/00_manuscript/state_vs_genotype_split/
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.celldrug import decompose, load_prism
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "00_manuscript"
TAB = ROOT / "results" / "tables"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"
ALPHAS = (100.0, 1000.0, 10000.0)


def cv_r2(X, y, folds, alphas=ALPHAS):
    """Dual ridge, alpha chosen inside the training folds (n < p everywhere)."""
    y = np.asarray(y, float)
    ss = float(np.sum((y - y.mean()) ** 2))
    if ss <= 0 or X.shape[1] == 0:
        return np.nan
    pred = np.zeros_like(y)
    for f in np.unique(folds):
        tr, te = folds != f, folds == f
        A, B = X[tr], X[te]
        mu, sd = A.mean(0), A.std(0)
        sd = np.where(sd > 0, sd, 1.0)
        A, B = (A - mu) / sd, (B - mu) / sd
        ytr = y[tr]
        inner = folds[tr] % 3
        best, best_s = alphas[0], -np.inf
        for al in alphas:
            sc = 0.0
            for g in np.unique(inner):
                it, iv = inner != g, inner == g
                Ka = A[it] @ A[it].T + al * np.eye(int(it.sum()))
                m = ytr[it].mean()
                p = m + (A[iv] @ A[it].T) @ np.linalg.solve(Ka, ytr[it] - m)
                sc -= float(np.sum((ytr[iv] - p) ** 2))
            if sc > best_s:
                best_s, best = sc, al
        Ka = A @ A.T + best * np.eye(A.shape[0])
        m = ytr.mean()
        pred[te] = m + (B @ A.T) @ np.linalg.solve(Ka, ytr - m)
    return 1 - float(np.sum((y - pred) ** 2)) / ss


def cluster_ci(vals, groups, n_boot=2000, seed=0):
    """Bootstrap resampling whole compounds, not compound-line observations.

    Compound residuals correlate at about r = 0.23, so an observation-level
    interval is far too narrow; the cluster is the compound.
    """
    rng = np.random.default_rng(seed)
    vals = np.asarray(vals, float)
    g = np.asarray(groups)
    ug = np.unique(g)
    idx = {u: np.where(g == u)[0] for u in ug}
    out = []
    for _ in range(n_boot):
        pick = rng.choice(ug, len(ug), replace=True)
        take = np.concatenate([idx[u] for u in pick])
        out.append(np.nanmedian(vals[take]))
    return float(np.nanmedian(vals)), float(np.percentile(out, 2.5)), \
        float(np.percentile(out, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-compounds", type=int, default=150)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(ROOT / "scripts"))
    from expression_gap_closure import load_expression

    R, K, lines, ti = load_prism()
    print(f"{R.shape[0]} lines x {R.shape[1]} conditions", flush=True)
    dec = decompose(R, K, lines)
    print(f"{len(dec.gamma_rep)} compounds with replicate-validated gamma",
          flush=True)

    s = dec.alpha_split.dropna()
    r_alpha = float(np.corrcoef(s.half_A, s.half_B)[0, 1])
    print(f"\ngeneral sensitivity reproduces across disjoint compound halves: "
          f"r = {r_alpha:.3f} (n = {len(s)} lines)")

    EXPR = load_expression()
    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(r"DepMap;\s*(ACH-\d+)")[0]
    cv2dep = dict(zip(smp.Cellosaurus_ID, dep))
    mut = pd.read_csv(PR / "mutation_long.tsv.gz", sep="\t", low_memory=False)
    mut["depmap"] = mut.RRID.map(cv2dep)
    ns = mut[mut.depmap.notna() &
             mut.Variant_Classification.astype(str).str.contains(
                 "Missense|Nonsense|Frame_Shift|Splice|In_Frame", na=False)]
    gc = ns.groupby("Gene_symbol").depmap.nunique()
    genes = sorted(gc[gc >= 20].index)
    gsets = {g: set(d.depmap) for g, d in
             ns[ns.Gene_symbol.isin(genes)].groupby("Gene_symbol")}
    ci = pd.read_csv(PR / "secondary-screen-cell-line-info.csv")
    tis = dict(zip(ci.depmap_id, ci.primary_tissue.astype(str)))
    print(f"predictors: {EXPR.shape[1]} expression features, {len(genes)} "
          f"mutated genes", flush=True)

    def blocks(keep):
        XE = EXPR.loc[keep].to_numpy(float)
        XE = XE[:, XE.std(0) > 0]
        XN = np.stack([np.isin(keep, list(gsets[g])).astype(float)
                       for g in genes], 1)
        XN = XN[:, XN.std(0) > 0]
        lin = np.array([tis.get(d, "?") for d in keep])
        ul = [u for u in np.unique(lin) if (lin == u).sum() >= 5]
        XL = (np.stack([(lin == u).astype(float) for u in ul], 1) if ul
              else np.zeros((len(keep), 0)))
        return {"expression": XE, "mutations": XN, "lineage": XL}

    # ---------- A. the cell property ----------
    keep_a = [d for d in s.index if d in EXPR.index]
    ya = s.loc[keep_a].mean(axis=1).to_numpy()
    ya = ya - ya.mean()
    folds_a = np.random.default_rng(0).permutation(np.arange(len(keep_a)) % 5)
    Ba = blocks(keep_a)
    print(f"\nA. predicting GENERAL SENSITIVITY (a cell property), "
          f"n = {len(keep_a)} lines")
    prop = {}
    for lab, X in Ba.items():
        prop[lab] = cv_r2(X, ya, folds_a)
        print(f"   {lab:12s} CV R2 {prop[lab]:+.4f}")

    # ---------- B. the cell-drug relation ----------
    print(f"\nB. predicting the CELL-DRUG RELATION, per compound")
    sel = sorted(dec.gamma, key=lambda c: -len(dec.gamma[c]))[:args.n_compounds]
    rows = []
    for i, cpd in enumerate(sel):
        y0 = dec.gamma[cpd]
        keep = [d for d in y0.index if d in EXPR.index]
        if len(keep) < 80:
            continue
        yv = y0.loc[keep].to_numpy(float)
        yv = yv - yv.mean()
        folds = np.random.default_rng(i).permutation(np.arange(len(keep)) % 5)
        Bk = blocks(keep)
        rec = {"compound": cpd, "n_lines": len(keep)}
        for lab, X in Bk.items():
            rec[lab] = cv_r2(X, yv, folds)
        rows.append(rec)
        if i % 25 == 0:
            print(f"   {i+1}/{len(sel)} {cpd[:20]:20s} expr "
                  f"{rec['expression']:+.4f}  mut {rec['mutations']:+.4f}",
                  flush=True)
    B = pd.DataFrame(rows)
    B.to_csv(TAB / "state_vs_genotype_split.csv", index=False)

    print(f"\n   compound-cluster bootstrap over {len(B)} compounds:")
    rel = {}
    for lab in ("expression", "mutations", "lineage"):
        m, lo, hi = cluster_ci(B[lab].to_numpy(), B.compound.to_numpy())
        rel[lab] = (m, lo, hi)
        flag = "" if (lo > 0) else "   <- interval includes zero"
        print(f"   {lab:12s} median CV R2 {m:+.4f}  95% CI "
              f"[{lo:+.4f}, {hi:+.4f}]{flag}")

    print("\nWhat this settles")
    pe, re_ = prop["expression"], rel["expression"][0]
    print(f"  expression predicts the cell property at R2 = {pe:+.4f}")
    print(f"  expression predicts the cell-drug relation at R2 = {re_:+.4f}"
          f"  ({re_/pe:.0%} of the property figure)" if pe > 0 else "")
    if rel["expression"][1] <= 0:
        print("  The relation term is NOT predicted by expression once general\n"
              "  sensitivity is removed. The published claim conflated the two.")
    else:
        print("  The relation term survives with expression above zero, so the\n"
              "  claim stands on a residual that is no longer general "
              "sensitivity.")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    ax[0].scatter(s.half_A, s.half_B, s=7, alpha=0.35, color=VIOLET,
                  edgecolors="none")
    ax[0].set_xlabel("general sensitivity, compound half A")
    ax[0].set_ylabel("general sensitivity, compound half B")
    ax[0].set_title(f"a  General sensitivity is a real cell property\n"
                    f"r = {r_alpha:.3f}, n = {len(s)} lines", loc="left",
                    fontweight="bold", fontsize=9.5)

    labs = ["expression", "lineage", "mutations"]
    xx = np.arange(3)
    ax[1].bar(xx - 0.19, [prop[l] for l in labs], width=0.36, color=GREY,
              label="cell property (general sensitivity)")
    ax[1].bar(xx + 0.19, [rel[l][0] for l in labs], width=0.36, color=ORANGE,
              label="cell-drug relation (interaction)")
    ax[1].errorbar(xx + 0.19, [rel[l][0] for l in labs],
                   yerr=[[rel[l][0] - rel[l][1] for l in labs],
                         [rel[l][2] - rel[l][0] for l in labs]],
                   fmt="none", ecolor="#333", lw=1.1, capsize=2.5)
    ax[1].axhline(0, color="#444", lw=0.9)
    ax[1].set_xticks(xx, labs, fontsize=8)
    ax[1].set_ylabel("cross-validated $R^2$")
    ax[1].legend(frameon=False, fontsize=7.2)
    ax[1].set_title("b  Which of the two does each block predict?", loc="left",
                    fontweight="bold", fontsize=9.5)

    ax[2].hist(B.expression.dropna(), bins=28, alpha=0.8, color=VIOLET,
               label="expression")
    ax[2].hist(B.mutations.dropna(), bins=28, alpha=0.65, color=ORANGE,
               label="mutations")
    ax[2].axvline(0, color="#444", lw=1.0)
    ax[2].set_xlabel("CV $R^2$ for the interaction, per compound")
    ax[2].set_ylabel("compounds")
    ax[2].legend(frameon=False, fontsize=7.5)
    ax[2].set_title(f"c  Per compound, after removing\n    general sensitivity "
                    f"(n = {len(B)})", loc="left", fontweight="bold",
                    fontsize=9.5)
    fig.suptitle("Separating the cell property from the cell-drug relation",
                 fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "state_vs_genotype_split", FIG,
                    source_data={"per_compound": B,
                                 "general_sensitivity": s.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
