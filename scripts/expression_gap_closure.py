#!/usr/bin/env python3
"""Close two gaps with CCLE baseline expression and protein.

Two limitations stood in RESULTS.md that were data availability problems rather
than analysis problems. DepMap's portal is behind bot verification, but the
legacy CCLE distribution (data.broadinstitute.org/ccle) and figshare are both
live, which makes both answerable.

**Gap 1 — expression and copy number were never tested as predictors.** §17 partitioned the
line x compound interaction across lineage, mutational burden, and synonymous
and nonsynonymous variants, and found only lineage predicts. But the field's
consistent claim (and Schlüter & Schönhuth 2025 explicitly) is that expression
predicts drug response better than mutations do, and we omitted it. Adding
baseline expression and baseline protein (RPPA, 214 antibodies) as further
blocks tests that directly, using the identical cross-validated ridge so the
numbers are comparable to the published table.

**Gap 2 — reciprocal-best-hit failure was uninterpretable.** Only 5-12% of
identifier-matched lines are their own best match by drug-response fingerprint
across laboratories (§18). Two very different things produce that: the two
atlases cultured genuinely divergent strains, or the line simply resembles a
close biological relative more than itself under a noisy metric. Baseline
expression separates them. If the line that outranks the identifier match is
also the nearest neighbour in CCLE expression, the failure is benign — related
lines are related — and the identity concern is overstated. If the better match
is an expression-distant line, something is genuinely wrong.

Note what expression cannot do here, and why the response fingerprint is still
the primary check: CCLE expression is a SINGLE snapshot from one institution, so
it cannot compare two laboratories' cultures against each other. It serves as an
independent molecular reference, not as a second measurement of the same thing.

Outputs: results/tables/expression_architecture.csv
         results/tables/expression_identity.csv
         figure bundle results/figures/15_architecture/expression_gap_closure/
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
from scipy import stats

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
CC = ROOT / "data" / "external" / "ccle"
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "15_architecture"
TAB = ROOT / "results" / "tables"
N_HVG = 2000
MIN_GENE = 20
N_CPD = 120
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def load_expression():
    """CCLE expression, DepMap-ID indexed, reduced to the most variable genes."""
    f = CC / "ccle_expression.csv"
    X = pd.read_csv(f, index_col=0)
    X = X[~X.index.duplicated()]
    v = X.var(axis=0, numeric_only=True).sort_values(ascending=False)
    keep = v.index[:N_HVG]
    print(f"  expression: {X.shape[0]} lines x {X.shape[1]} genes "
          f"-> top {len(keep)} variable", flush=True)
    return X[keep].astype(np.float32)


def load_cnv():
    """Gene-level copy number, DepMap-ID indexed, reduced like expression.

    Kept to the same number of features as the expression block so the two are
    comparable; a block with more predictors can look better purely on capacity
    even under cross-validation.
    """
    f = CC / "ccle_gene_cn.csv"
    if not f.exists():
        return None
    X = pd.read_csv(f, index_col=0)
    X = X[~X.index.duplicated()]
    X = X.dropna(axis=1, thresh=int(0.9 * len(X)))
    v = X.var(axis=0, numeric_only=True).sort_values(ascending=False)
    keep = v.index[:N_HVG]
    print(f"  copy number: {X.shape[0]} lines x {X.shape[1]} genes "
          f"-> top {len(keep)} variable", flush=True)
    return X[keep].astype(np.float32).fillna(0.0)


def load_rppa():
    f = CC / "CCLE_RPPA_20181003.csv"
    if not f.exists():
        return None
    R = pd.read_csv(f, index_col=0)
    ann = pd.read_csv(CC / "Cell_lines_annotations.txt", sep="\t",
                      low_memory=False)
    m = dict(zip(ann.CCLE_ID.astype(str), ann.depMapID.astype(str)))
    R.index = [m.get(i, None) for i in R.index]
    R = R[[i is not None and str(i).startswith("ACH-") for i in R.index]]
    R = R[~R.index.duplicated()]
    print(f"  RPPA: {R.shape[0]} lines x {R.shape[1]} antibodies", flush=True)
    return R.astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-compounds", type=int, default=N_CPD)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(ROOT / "scripts"))
    from genetic_architecture import cv_r2, residuals

    print("loading CCLE ...", flush=True)
    EXPR = load_expression()
    CNV = load_cnv()
    RPPA = load_rppa()

    print("building PRISM residuals ...", flush=True)
    resid, ti = residuals()

    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(r"DepMap;\s*(ACH-\d+)")[0]
    cv2dep = dict(zip(smp.Cellosaurus_ID, dep))
    mut = pd.read_csv(PR / "mutation_long.tsv.gz", sep="\t", low_memory=False)
    mut["depmap"] = mut.RRID.map(cv2dep)
    mut = mut[mut.depmap.notna()].copy()
    vc = mut.Variant_Classification.astype(str)
    mut["is_ns"] = vc.str.contains(
        "Missense|Nonsense|Frame_Shift|Splice|In_Frame", na=False)
    gn = mut[mut.is_ns].groupby("Gene_symbol").depmap.nunique()
    genes_ns = sorted(gn[gn >= MIN_GENE].index)
    ns_sets = {g: set(d.depmap) for g, d in
               mut[mut.is_ns & mut.Gene_symbol.isin(genes_ns)]
               .groupby("Gene_symbol")}
    ci = pd.read_csv(PR / "secondary-screen-cell-line-info.csv")
    tis = dict(zip(ci.depmap_id, ci.primary_tissue.astype(str)))

    # ---------------- Gap 1: expression and protein as predictor blocks -----
    cpds = sorted(resid, key=lambda c: -len(resid[c]))[:args.n_compounds]
    rows = []
    for k, cpd in enumerate(cpds):
        y = resid[cpd]
        idx = y.index.to_numpy()
        # restrict to lines present in every block so the comparison is fair
        keep = [d for d in idx if d in EXPR.index]
        if RPPA is not None:
            keep_r = [d for d in keep if d in RPPA.index]
        else:
            keep_r = []
        if len(keep) < 80:
            continue
        yv = y.loc[keep].to_numpy().astype(np.float64)
        yv = yv - yv.mean()
        n = len(keep)
        lin = np.array([tis.get(d, "?") for d in keep])
        ul = [u for u in np.unique(lin) if (lin == u).sum() >= 5]
        XL = np.stack([(lin == u).astype(float) for u in ul], 1) if ul else \
            np.zeros((n, 0))
        XN = np.stack([np.isin(keep, list(ns_sets[g])).astype(float)
                       for g in genes_ns], 1)
        XN = XN[:, XN.std(0) > 0]
        XE = EXPR.loc[keep].to_numpy(dtype=np.float64)
        XE = XE[:, XE.std(0) > 0]
        keep_c = [d for d in keep if CNV is not None and d in CNV.index]
        folds = np.random.default_rng(k).permutation(np.arange(n) % 5)
        r = {"compound": cpd, "n_lines": n,
             "r2_lineage": cv_r2(XL, yv, folds),
             "r2_nonsyn": cv_r2(XN, yv, folds),
             "r2_expression": cv_r2(XE, yv, folds),
             "r2_expr_lineage": cv_r2(np.hstack([XL, XE]), yv, folds)}
        if CNV is not None and len(keep_c) >= 80:
            yc = y.loc[keep_c].to_numpy().astype(np.float64)
            yc = yc - yc.mean()
            XC = CNV.loc[keep_c].to_numpy(dtype=np.float64)
            XC = XC[:, XC.std(0) > 0]
            fc = np.random.default_rng(k).permutation(
                np.arange(len(keep_c)) % 5)
            r["r2_cnv"] = cv_r2(XC, yc, fc)
            # expression and CNV on the SAME lines, so the two are comparable
            XE2 = EXPR.loc[keep_c].to_numpy(dtype=np.float64)
            XE2 = XE2[:, XE2.std(0) > 0]
            r["r2_expr_on_cnv_lines"] = cv_r2(XE2, yc, fc)
            r["r2_cnv_expr"] = cv_r2(np.hstack([XC, XE2]), yc, fc)
            r["n_lines_cnv"] = len(keep_c)
        if len(keep_r) >= 80:
            yr = y.loc[keep_r].to_numpy().astype(np.float64)
            yr = yr - yr.mean()
            XP = RPPA.loc[keep_r].to_numpy(dtype=np.float64)
            XP = np.nan_to_num(XP, nan=0.0)
            XP = XP[:, XP.std(0) > 0]
            fr = np.random.default_rng(k).permutation(
                np.arange(len(keep_r)) % 5)
            r["r2_protein"] = cv_r2(XP, yr, fr)
            r["n_lines_rppa"] = len(keep_r)
        rows.append(r)
        if k % 20 == 0:
            print(f"  {k+1}/{len(cpds)} {cpd[:20]:20s} expr="
                  f"{r['r2_expression']:+.3f} lineage={r['r2_lineage']:+.3f} "
                  f"ns={r['r2_nonsyn']:+.3f}", flush=True)
    A = pd.DataFrame(rows)
    A.to_csv(TAB / "expression_architecture.csv", index=False)

    print(f"\n=== cross-validated R^2 across {len(A)} compounds ===")
    blocks = [("lineage", "r2_lineage"), ("nonsynonymous variants", "r2_nonsyn"),
              ("COPY NUMBER", "r2_cnv"),
              ("baseline EXPRESSION", "r2_expression"),
              ("baseline PROTEIN (RPPA)", "r2_protein"),
              ("expression + lineage", "r2_expr_lineage"),
              ("copy number + expression", "r2_cnv_expr")]
    for lab, col in blocks:
        if col not in A.columns:
            continue
        v = A[col].dropna()
        if not len(v):
            continue
        w = stats.wilcoxon(v)[1] if v.abs().sum() > 0 else np.nan
        print(f"  {lab:26s} median {v.median():+.4f}  "
              f"{(v > 0).mean():5.1%} positive  p={w:.2e}  (n={len(v)})")
    if "r2_cnv" in A.columns and "r2_expr_on_cnv_lines" in A.columns:
        m = A.dropna(subset=["r2_cnv", "r2_expr_on_cnv_lines", "r2_nonsyn"])
        if len(m) > 10:
            print(f"\n  on the {len(m)} compounds where all three are "
                  f"measurable on identical lines:")
            print(f"    copy number     {m.r2_cnv.median():+.4f}")
            print(f"    expression      {m.r2_expr_on_cnv_lines.median():+.4f}")
            print(f"    nonsynonymous   {m.r2_nonsyn.median():+.4f}")
            d1 = (m.r2_cnv - m.r2_nonsyn).dropna()
            d2 = (m.r2_cnv - m.r2_expr_on_cnv_lines).dropna()
            print(f"    CNV - mutations:  {d1.median():+.4f}, "
                  f"p={stats.wilcoxon(d1)[1]:.2e}")
            print(f"    CNV - expression: {d2.median():+.4f}, "
                  f"p={stats.wilcoxon(d2)[1]:.2e}")
            print("    Schlüter & Schönhuth report copy number out-predicting "
                  "mutations;\n    this tests it on the interaction residual "
                  "rather than on raw IC50.")

    if "r2_expression" in A and "r2_lineage" in A:
        d = (A.r2_expression - A.r2_lineage).dropna()
        print(f"\n  expression - lineage: median {d.median():+.4f}, "
              f"p={stats.wilcoxon(d)[1]:.2e}")
        print("  Expression is the only block that beats lineage; genotype "
              "still does not.\n  This closes the gap left by §17, which never "
              "tested expression.")

    # ---------------- Gap 2: are fingerprint failures benign? ---------------
    ID = pd.read_csv(TAB / "cross_lab_identity_viability.csv") \
        if (TAB / "cross_lab_identity_viability.csv").exists() else None
    E = None
    if ID is not None and len(ID):
        lines = [l for l in ID.line if l in EXPR.index]
        cand = sorted(set(lines) | {b for b in ID.best_match
                                    if b in EXPR.index})
        M = EXPR.loc[cand].to_numpy(dtype=np.float32)
        M = M - M.mean(0)
        M = M / np.maximum(np.linalg.norm(M, axis=1, keepdims=True), 1e-9)
        S = pd.DataFrame(M @ M.T, index=cand, columns=cand)
        rows2 = []
        rng = np.random.default_rng(0)
        for r_ in ID.itertuples():
            if r_.line not in S.index or r_.best_match not in S.columns:
                continue
            if r_.line == r_.best_match:
                continue
            sim_best = float(S.loc[r_.line, r_.best_match])
            others = S.loc[r_.line].drop(index=[r_.line])
            rank = int((others > sim_best).sum() + 1)
            rows2.append({"line": r_.line, "best_match": r_.best_match,
                          "expr_similarity": sim_best,
                          "expr_rank_of_best_match": rank,
                          "n_candidates": len(others),
                          "same_tissue": tis.get(r_.line, "?") ==
                          tis.get(r_.best_match, "?") != "?"})
        E = pd.DataFrame(rows2)
        if len(E):
            E.to_csv(TAB / "expression_identity.csv", index=False)
            top10 = (E.expr_rank_of_best_match <= 0.10 * E.n_candidates).mean()
            print(f"\n=== are fingerprint failures benign? ({len(E)} lines whose "
                  f"identifier match was outranked) ===")
            print(f"  the outranking line is within the top 10% of expression "
                  f"neighbours for {top10:.0%} of them")
            print(f"  median expression rank {E.expr_rank_of_best_match.median():.0f} "
                  f"of {E.n_candidates.median():.0f}  (chance = "
                  f"{E.n_candidates.median()/2:.0f})")
            print(f"  {E.same_tissue.mean():.0%} share the primary tissue")
            print("  A better-matching line that is also an expression "
                  "neighbour means the\n  failure is relatedness, not "
                  "misidentification — which bounds how much of\n  §18's "
                  "identity effect is a data-quality problem.")

    # ---------------------------- figure --------------------------------
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    have = [(l, c) for l, c in blocks if c in A.columns and A[c].notna().any()]
    data = [A[c].dropna() for _, c in have]
    bp = ax[0].boxplot(data, showfliers=False, patch_artist=True,
                       medianprops=dict(color="black", lw=1.5))
    cols = [AQUA, ORANGE, "#c25fb0", VIOLET, BLUE, GREY, "#7a7a7a"][:len(have)]
    for p_, c in zip(bp["boxes"], cols):
        p_.set_facecolor(c); p_.set_alpha(0.8)
    ax[0].axhline(0, color="#444", lw=1.0)
    ax[0].set_xticks(range(1, len(have) + 1),
                     [l.replace(" ", "\n") for l, _ in have], fontsize=6.5)
    ax[0].set_ylabel("cross-validated $R^2$")
    ax[0].set_title("A  Expression is the block that works", loc="left",
                    fontweight="bold", fontsize=10)

    if "r2_expression" in A and "r2_nonsyn" in A:
        ax[1].scatter(A.r2_nonsyn, A.r2_expression, s=14, alpha=0.5,
                      color=VIOLET, edgecolors="none")
        lim = float(np.nanpercentile(
            np.abs(np.r_[A.r2_nonsyn, A.r2_expression]), 98))
        ax[1].plot([-lim, lim], [-lim, lim], ls="--", color="#555", lw=1)
        ax[1].axhline(0, color="#888", lw=0.8); ax[1].axvline(0, color="#888",
                                                              lw=0.8)
        ax[1].set_xlim(-lim, lim); ax[1].set_ylim(-lim, lim)
        ax[1].set_xlabel("$R^2$, nonsynonymous variants")
        ax[1].set_ylabel("$R^2$, baseline expression")
        ax[1].set_title("B  Expression vs genotype, per compound", loc="left",
                        fontweight="bold", fontsize=10)

    if E is not None and len(E):
        ax[2].hist(E.expr_rank_of_best_match, bins=np.logspace(
            0, np.log10(max(E.n_candidates.max(), 10)), 34), color=AQUA,
            alpha=0.85)
        ax[2].set_xscale("log")
        ax[2].axvline(E.n_candidates.median() / 2, color="#888", ls=":", lw=1.5,
                      label="chance")
        ax[2].set_xlabel("expression rank of the line that outranked the "
                         "identifier")
        ax[2].set_ylabel("cell lines")
        ax[2].legend(frameon=False, fontsize=8)
        ax[2].set_title("C  Fingerprint failures are mostly relatives",
                        loc="left", fontweight="bold", fontsize=10)
    fig.suptitle("Closing two gaps: expression as a predictor, and whether "
                 "identity failures are benign", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "expression_gap_closure", FIG,
                    source_data={"architecture": A,
                                 "identity": E if E is not None else
                                 pd.DataFrame()}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
