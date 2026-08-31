#!/usr/bin/env python3
"""The genetic architecture of drug-response variation between cellular contexts.

This reframes the project's question into the one Nature Genetics actually asks.
So far we have asked "can genotype predict context-dependent drug response",
found the known biomarkers, and stopped. The architecture question is different
and has not been answered for this phenotype: **how much of the variation between
cellular contexts is genetic at all, is it oligogenic or polygenic, and how much
of what looks genetic is population structure rather than mechanism?**

The design rests on one control that makes the whole thing interpretable.

**Synonymous variants are the ideal negative control.** A synonymous change
cannot alter the protein, so it cannot mechanistically cause a difference in drug
response. But it carries exactly the same population structure, ancestry and
lineage information as a nonsynonymous variant in the same gene, and — because
CCLE calls variants without matched normals — the same germline contamination.
So the comparison

    predictive power of NONSYNONYMOUS variants   (mechanism + confounding)
    predictive power of SYNONYMOUS variants      (confounding only)

isolates the mechanistic component. If the two are equal, apparent "genotype
predicts drug response" is population structure. If nonsynonymous clearly
exceeds synonymous, the excess is a genuine, quantified mechanistic signal. This
is the drug-response analogue of using synonymous sites as a neutral reference in
molecular evolution, and it is what the §15 germline finding demands: three of
the four candidate alleles sat beside recurrent SYNONYMOUS variants, which is how
we caught them.

Four predictor blocks per cell line:

  lineage       one-hot primary tissue — the non-genetic context baseline
  synonymous    silent variants in frequently mutated genes — confounding only
  nonsynonymous missense/nonsense/frameshift in the same genes — mechanism + confounding
  burden        total mutation count — a single summary covariate

Everything is scored by **cross-validated** R² (5-fold, ridge), because in-sample
R² rises with predictor count and would make the larger block win by construction.
Blocks are also fitted jointly to report what each adds over lineage, since
genotype and lineage are strongly correlated and the marginal contribution is the
quantity of interest.

Outputs: results/tables/genetic_architecture.csv
         results/tables/genetic_architecture_summary.csv
         figure bundle results/figures/15_architecture/genetic_architecture/
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
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "15_architecture"
TAB = ROOT / "results" / "tables"
MIN_LINES = 100
MIN_GENE = 20
N_CPD = 300
ALPHAS = np.logspace(-1, 4, 12)
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def cv_r2(Xb, y, folds, alphas=ALPHAS):
    """Out-of-sample R^2 of ridge on a predictor block, alpha chosen inside
    each training fold so the estimate stays honest."""
    if Xb.shape[1] == 0:
        return 0.0
    pred = np.zeros_like(y)
    for f in np.unique(folds):
        tr, te = folds != f, folds == f
        Xtr, Xte = Xb[tr], Xb[te]
        mu, sd = Xtr.mean(0), Xtr.std(0)
        sd = np.where(sd > 0, sd, 1.0)
        Xtr = (Xtr - mu) / sd
        Xte = (Xte - mu) / sd
        ytr = y[tr]
        ym = ytr.mean()

        def fit_predict(A, ya, B, a):
            """Ridge in whichever space is smaller.

            With ~740 cell lines and up to ~12,000 gene predictors the primal
            p x p solve is hopeless; the dual form solves an n x n system and
            gives identical predictions, so the block size stops driving the
            cost.
            """
            if A.shape[1] <= A.shape[0]:
                w = np.linalg.solve(A.T @ A + a * np.eye(A.shape[1]), A.T @ ya)
                return B @ w
            Kk = A @ A.T
            al = np.linalg.solve(Kk + a * np.eye(A.shape[0]), ya)
            return (B @ A.T) @ al

        inner = np.arange(len(ytr)) % 3
        best, best_e = alphas[0], np.inf
        for a in alphas:
            err = 0.0
            for g in range(3):
                itr, ite = inner != g, inner == g
                im = ytr[itr].mean()
                pr_ = fit_predict(Xtr[itr], ytr[itr] - im, Xtr[ite], a)
                err += float(np.sum((ytr[ite] - im - pr_) ** 2))
            if err < best_e:
                best_e, best = err, a
        pred[te] = ym + fit_predict(Xtr, ytr - ym, Xte, best)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def residuals():
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
    ap.add_argument("--n-compounds", type=int, default=N_CPD)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    resid, ti = residuals()
    print(f"{len(resid)} compounds with interaction residuals", flush=True)

    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    import re
    dep = smp.cross_references.astype(str).str.extract(r"DepMap;\s*(ACH-\d+)")[0]
    cv2dep = dict(zip(smp.Cellosaurus_ID, dep))
    mut = pd.read_csv(PR / "mutation_long.tsv.gz", sep="\t", low_memory=False)
    mut["depmap"] = mut.RRID.map(cv2dep)
    mut = mut[mut.depmap.notna()].copy()
    vc = mut.Variant_Classification.astype(str)
    mut["is_syn"] = vc.str.contains("Silent|Synonymous", na=False)
    mut["is_ns"] = vc.str.contains(
        "Missense|Nonsense|Frame_Shift|Splice|In_Frame", na=False)
    burden = mut.groupby("depmap").size()

    # genes carrying BOTH classes often enough to be compared fairly
    gs = mut[mut.is_syn].groupby("Gene_symbol").depmap.nunique()
    gn = mut[mut.is_ns].groupby("Gene_symbol").depmap.nunique()
    # The two blocks must be matched or the comparison is meaningless: an
    # 11,809-gene nonsynonymous block would beat a 3,447-gene synonymous one on
    # predictor count alone. Restrict BOTH to the genes that carry each class in
    # >=MIN_GENE lines, so the blocks cover the same genes and are the same size,
    # and the only difference is whether the variant can change the protein.
    both = sorted(set(gs[gs >= MIN_GENE].index) & set(gn[gn >= MIN_GENE].index))
    genes_syn = genes_ns = both
    print(f"{(gn >= MIN_GENE).sum()} genes nonsynonymous-mutated and "
          f"{(gs >= MIN_GENE).sum()} synonymous-mutated in >={MIN_GENE} lines; "
          f"using the {len(both)} carrying BOTH so the blocks are matched",
          flush=True)

    ci = pd.read_csv(PR / "secondary-screen-cell-line-info.csv")
    tis = dict(zip(ci.depmap_id, ci.primary_tissue.astype(str)))

    syn_sets = {g: set(d.depmap) for g, d in
                mut[mut.is_syn & mut.Gene_symbol.isin(genes_syn)]
                .groupby("Gene_symbol")}
    ns_sets = {g: set(d.depmap) for g, d in
               mut[mut.is_ns & mut.Gene_symbol.isin(genes_ns)]
               .groupby("Gene_symbol")}

    cpds = sorted(resid, key=lambda c: -len(resid[c]))[:args.n_compounds]
    rows = []
    for k, cpd in enumerate(cpds):
        y = resid[cpd]
        idx = y.index.to_numpy()
        yv = y.to_numpy().astype(np.float64)
        yv = yv - yv.mean()
        n = len(idx)
        lin = np.array([tis.get(d, "?") for d in idx])
        ul = [u for u in np.unique(lin) if (lin == u).sum() >= 5]
        XL = np.stack([(lin == u).astype(float) for u in ul], 1) if ul else \
            np.zeros((n, 0))
        XS = np.stack([np.isin(idx, list(syn_sets[g])).astype(float)
                       for g in genes_syn], 1) if genes_syn else np.zeros((n, 0))
        XN = np.stack([np.isin(idx, list(ns_sets[g])).astype(float)
                       for g in genes_ns], 1) if genes_ns else np.zeros((n, 0))
        XB = burden.reindex(idx).fillna(0).to_numpy().reshape(-1, 1).astype(float)
        XB = np.log1p(XB)
        # drop invariant columns within this compound's line set
        XS = XS[:, XS.std(0) > 0]
        XN = XN[:, XN.std(0) > 0]
        folds = np.arange(n) % 5
        rng = np.random.default_rng(k)
        folds = rng.permutation(folds)
        r = {"compound": cpd, "n_lines": n,
             "n_syn": XS.shape[1], "n_ns": XN.shape[1]}
        r["r2_lineage"] = cv_r2(XL, yv, folds)
        r["r2_burden"] = cv_r2(XB, yv, folds)
        r["r2_syn"] = cv_r2(XS, yv, folds)
        r["r2_ns"] = cv_r2(XN, yv, folds)
        r["r2_lin_syn"] = cv_r2(np.hstack([XL, XS]), yv, folds)
        r["r2_lin_ns"] = cv_r2(np.hstack([XL, XN]), yv, folds)
        r["r2_all"] = cv_r2(np.hstack([XL, XS, XN, XB]), yv, folds)
        rows.append(r)
        if k % 25 == 0:
            print(f"  {k+1}/{len(cpds)} {cpd[:22]:22s} lineage="
                  f"{r['r2_lineage']:+.3f} syn={r['r2_syn']:+.3f} "
                  f"ns={r['r2_ns']:+.3f}", flush=True)
    A = pd.DataFrame(rows)
    A["ns_minus_syn"] = A.r2_ns - A.r2_syn
    A["ns_over_lineage"] = A.r2_lin_ns - A.r2_lineage
    A["syn_over_lineage"] = A.r2_lin_syn - A.r2_lineage
    A["mechanistic"] = A.ns_over_lineage - A.syn_over_lineage
    moa = ti.drop_duplicates("name").set_index("name").moa
    A["moa"] = A.compound.map(moa)
    A.to_csv(TAB / "genetic_architecture.csv", index=False)

    def line(lbl, col):
        v = A[col]
        w = stats.wilcoxon(v)[1] if v.abs().sum() > 0 else np.nan
        print(f"  {lbl:44s} median {v.median():+.4f}  "
              f"{(v > 0).mean():5.1%} positive  p={w:.2e}")

    print(f"\n=== cross-validated R^2 across {len(A)} compounds ===")
    for lbl, col in (("lineage alone", "r2_lineage"),
                     ("mutational burden alone", "r2_burden"),
                     ("SYNONYMOUS variants alone (confounding only)", "r2_syn"),
                     ("NONSYNONYMOUS variants alone (mech + confounding)",
                      "r2_ns"),
                     ("everything", "r2_all")):
        line(lbl, col)
    print("\n=== the mechanistic excess ===")
    line("nonsynonymous - synonymous (marginal)", "ns_minus_syn")
    line("nonsynonymous added over lineage", "ns_over_lineage")
    line("synonymous added over lineage", "syn_over_lineage")
    line("MECHANISTIC = (ns|lineage) - (syn|lineage)", "mechanistic")
    w = stats.wilcoxon(A.r2_ns, A.r2_syn)
    print(f"\n  paired nonsynonymous vs synonymous: p={w.pvalue:.2e}")
    print("  If these are indistinguishable, the apparent genetic signal in "
          "drug response\n  is population structure, not mechanism — "
          "synonymous variants cannot cause it.")

    S = pd.DataFrame({
        "block": ["lineage", "burden", "synonymous", "nonsynonymous", "all"],
        "median_cv_r2": [A.r2_lineage.median(), A.r2_burden.median(),
                         A.r2_syn.median(), A.r2_ns.median(), A.r2_all.median()],
        "frac_positive": [(A.r2_lineage > 0).mean(), (A.r2_burden > 0).mean(),
                          (A.r2_syn > 0).mean(), (A.r2_ns > 0).mean(),
                          (A.r2_all > 0).mean()]})
    S.to_csv(TAB / "genetic_architecture_summary.csv", index=False)
    print("\n" + S.round(4).to_string(index=False))

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    blocks = ["r2_lineage", "r2_burden", "r2_syn", "r2_ns", "r2_all"]
    labs = ["lineage", "burden", "synonymous\n(control)", "nonsynonymous",
            "all"]
    ax[0].boxplot([A[b].clip(-0.2, 0.6) for b in blocks], showfliers=False,
                  medianprops=dict(color=ORANGE, lw=2))
    ax[0].set_xticks(range(1, 6), labs, fontsize=7.5)
    ax[0].axhline(0, color="#444444", lw=0.9)
    ax[0].set_ylabel("cross-validated R² of the interaction residual")
    ax[0].set_title("A  What predicts context-dependent response",
                    loc="left", fontweight="bold", fontsize=10)

    lim = float(np.nanpercentile(np.abs(np.r_[A.r2_syn, A.r2_ns]), 98))
    ax[1].scatter(A.r2_syn, A.r2_ns, s=10, alpha=0.4, color=BLUE,
                  edgecolors="none")
    ax[1].plot([-lim, lim], [-lim, lim], ls="--", color="#555555", lw=1)
    ax[1].set_xlim(-lim, lim); ax[1].set_ylim(-lim, lim)
    ax[1].set_xlabel("R², synonymous (confounding only)")
    ax[1].set_ylabel("R², nonsynonymous")
    ax[1].set_title(f"B  Above the line = mechanism\n"
                    f"median excess {A.ns_minus_syn.median():+.4f}",
                    loc="left", fontweight="bold", fontsize=10)

    ax[2].hist(A.mechanistic.clip(-0.15, 0.15), bins=40, color=VIOLET,
               alpha=0.85)
    ax[2].axvline(0, color="#444444", lw=1.2)
    ax[2].axvline(A.mechanistic.median(), color=ORANGE, lw=2)
    ax[2].set_xlabel("mechanistic component of R² (ns − syn, over lineage)")
    ax[2].set_ylabel("compounds")
    ax[2].set_title("C  How much is genuinely mechanistic", loc="left",
                    fontweight="bold", fontsize=10)
    fig.suptitle("Genetic architecture of drug-response variation between "
                 "cellular contexts", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "genetic_architecture", FIG,
                    source_data={"per_compound": A, "summary": S},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
