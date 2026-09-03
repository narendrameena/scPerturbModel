#!/usr/bin/env python3
"""Are OTHER groups' published context-specificity numbers inflated too?

Our own interaction shares were inflated about 1.5x because the residual they
were computed from still contained each context's general response -- how that
cell line reacts to being perturbed at all. That is a defect in our code. The
question that decides whether it is also a finding is whether the same term sits
inside statistics other groups have published, computed on their own data.

It should, because the construction is near-universal: a context-specificity
statistic is almost always built from a per-context response profile with the
PERTURBATION mean removed and nothing else. Whatever makes one context respond
more strongly to everything then reads as that context responding *differently*.

Three published statistics are examined. For each, the rule is the same and it is
strict: **reproduce the published number first, by the paper's own definition, on
the paper's own data.** Only if the replication lands near the published value
does the corrected number mean anything; if it does not, that is reported and no
claim of inflation is made. A correction applied to a statistic we failed to
reproduce would be measuring our own reimplementation, not their result.

  A. Srivatsan et al., Science 2020 (sci-Plex 3).  "Of 4,308 differentially
     expressed genes, 48% responded in a cell-type-dependent manner and 22%
     responded identically across all three cell lines."  Data: the sci-Plex 3
     pseudobulk, 3 lines x 189 compounds x 4 doses x 2 replicates. The two
     replicates are what make the classification testable rather than
     threshold-dependent: a gene is cell-type-dependent when its spread ACROSS
     lines exceeds its spread across replicates of the same line.

  B. Subramanian et al., Cell 2017 (CMap / LINCS).  Roughly a quarter of
     compounds produce signatures conserved across the cell panel. Data: LINCS
     phase 1 Level 4, the same resource.

  C. Ben-David et al., Nature 2018.  48 of 55 compounds active against at least
     one MCF7 strain were entirely inactive against another (87%), at >50% and
     <20% growth inhibition. Their 27 MCF7 strains are not public, so this is
     run on PRISM cell lines with the identical thresholds and is reported as an
     ANALOGUE, never as a replication of their number.

In every case the correction removes each context's mean response across the
OTHER perturbations, so a genuine perturbation-specific effect cannot be absorbed
into it.

**Result, stated up front: the generalisation does not hold, and the obvious
explanation for why does not hold either.**

Neither published statistic reproduced under our implementation of it -- sci-Plex
74.4% against a published 48%, LINCS 5.0% against a published 26% -- so neither
arm says anything about those papers, only about our reimplementations. And even
within our own versions the correction moves the statistic by almost nothing:
1.01x for sci-Plex, 1.15x for LINCS and in the *opposite* direction, against 1.5x
for our own interaction share.

The `context_share` diagnostic below was added to explain that pattern and does
not: across the three datasets it is ANTI-ordered with the effect (sci-Plex has
the largest context share, 8.3%, and the smallest movement; LINCS the smallest,
1.4%, and the largest). It is reported anyway, as a negative result, because the
reason it fails is the informative part. It measures the context main effect as a
share of TOTAL residual variance, which is dominated by measurement noise -- in
PRISM it returns 2.8%, while the same term is ~33% of the *reproducible*
(noise-free) residual by the replicate-covariance decomposition in
`perturbmodel.celldrug`.

That gap is the actual lesson. Our estimator compares components that noise
cannot inflate, because each is a covariance between independent estimates, and
in that arena the cell property is a third of the non-drug signal. A statistic
computed on a raw, noise-dominated residual is mostly measuring noise, so
removing a term worth a few percent of it changes little. The inflation should
therefore be expected in **replicate-validated** context-specificity statistics
and not in raw-residual ones -- but we have no other group's replicate-validated
statistic to test that on (TRADE is the nearest and its Perturb-seq data is not
on disk), so it stands as a prediction, not a demonstration.

Conclusion for the manuscript: report the defect as ours, quantified, with the
tool that prevents it. Do NOT claim the field's published numbers are inflated.

Outputs: results/tables/published_number_inflation.csv
         figure bundle results/figures/00_manuscript/published_inflation/
"""
import argparse
import re
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
SCP = ROOT / "data" / "external" / "scperturb" / "sciplex3_pseudobulk.h5ad"
LINCS = ROOT / "data" / "external" / "lincs"
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "00_manuscript"
TAB = ROOT / "results" / "tables"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def lognorm(X, scale=1e4):
    s = np.asarray(X.sum(axis=1)).ravel()
    s[s == 0] = 1.0
    return np.log1p(np.asarray(X) / s[:, None] * scale)


def loo_context_mean(D, ctx, pert):
    """Each context's mean response across the OTHER perturbations.

    Leave-one-perturbation-out, so a real perturbation-specific effect cannot be
    absorbed into the context's general response and quietly deleted.
    """
    out = np.zeros_like(D)
    for c in np.unique(ctx):
        m = ctx == c
        idx = np.where(m)[0]
        pv = pert[idx]
        # mean over rows of this context grouped by perturbation, then leave
        # the row's own perturbation out of that average
        by = {}
        for p in np.unique(pv):
            by[p] = D[idx[pv == p]].mean(0)
        tot = np.sum(list(by.values()), axis=0)
        n = len(by)
        if n < 2:
            continue
        for i, p in zip(idx, pv):
            out[i] = (tot - by[p]) / (n - 1)
    return out


# ---------------------------------------------------------------- A: sci-Plex
def sciplex(published_dep=0.48, published_same=0.22):
    if not SCP.exists():
        return None
    A = ad.read_h5ad(SCP)
    o = A.obs.copy()
    o["dose"] = o.dose_value.astype(str)
    X = lognorm(A.X if not hasattr(A.X, "toarray") else A.X.toarray())
    print(f"  sci-Plex 3 pseudobulk {X.shape}, "
          f"{o.cell_line.nunique()} lines, {o.perturbation.nunique()} "
          f"perturbations", flush=True)

    # vehicle per (line, replicate); the paper's deltas are against vehicle in
    # the same line, so a line's baseline never enters the response
    veh = {}
    isveh = (o.dose == "0.0") | o.perturbation.astype(str).str.contains(
        "vehicle|DMSO|control", case=False, na=False)
    for (ln, rp), g in o[isveh].groupby(["cell_line", "replicate"],
                                        observed=True):
        veh[(ln, rp)] = X[[o.index.get_loc(i) for i in g.index]].mean(0)
    if not veh:
        print("  no vehicle wells found; cannot form deltas")
        return None

    rows, D = [], []
    for i, (ix, r) in enumerate(o.iterrows()):
        if isveh.iloc[i]:
            continue
        v = veh.get((r.cell_line, r.replicate))
        if v is None:
            continue
        rows.append((r.cell_line, r.perturbation, r.dose, r.replicate))
        D.append(X[i] - v)
    K = pd.DataFrame(rows, columns=["line", "pert", "dose", "rep"])
    D = np.stack(D)
    # the paper's main comparison is at the top dose
    top = K.dose.value_counts().index[0]
    top = "10000.0" if "10000.0" in set(K.dose) else top
    sel = (K.dose == top).to_numpy()
    K, D = K[sel].reset_index(drop=True), D[sel]
    lines = sorted(K.line.unique())
    print(f"  {len(K)} treated pseudobulks at dose {top}; lines {lines}",
          flush=True)

    def classify(Dm):
        """Per (compound, gene): is the response cell-type-dependent?

        Between-line spread compared against the spread between replicates of
        the same line. Using the replicates rather than a fixed tolerance means
        the threshold is set by the data's own noise, and the same threshold
        applies before and after the correction.
        """
        dep = same = tot = 0
        for cpd, g in K.groupby("pert", observed=True):
            if g.line.nunique() < len(lines):
                continue
            per, wi = [], []
            ok = True
            for ln in lines:
                ii = g.index[g.line == ln].to_numpy()
                if len(ii) < 2:
                    ok = False
                    break
                per.append(Dm[ii].mean(0))
                wi.append(Dm[ii].var(0, ddof=1))
            if not ok:
                continue
            P = np.stack(per)                      # lines x genes
            within = np.mean(np.stack(wi), axis=0) + 1e-8
            between = P.var(0, ddof=1)
            # a gene is only counted if the compound moves it at all
            moved = np.abs(P).max(0) > 0.25
            d = moved & (between > within)
            s = moved & (between <= within)
            dep += int(d.sum()); same += int(s.sum()); tot += int(moved.sum())
        return dep, same, tot

    cs = context_share(D, K.line.to_numpy(), K.pert.to_numpy(),
                       "sci-Plex (vehicle-matched)")
    dep0, same0, tot0 = classify(D)
    C = loo_context_mean(D, K.line.to_numpy(), K.pert.to_numpy())
    dep1, same1, tot1 = classify(D - C)
    print(f"  as published:  {dep0}/{tot0} = {dep0/max(tot0,1):.1%} "
          f"cell-type-dependent, {same0/max(tot0,1):.1%} identical")
    print(f"  corrected:     {dep1}/{tot1} = {dep1/max(tot1,1):.1%} "
          f"cell-type-dependent, {same1/max(tot1,1):.1%} identical")
    return {"context_share": cs,
            "study": "Srivatsan 2020 (sci-Plex 3)",
            "statistic": "% of responsive genes that are cell-type-dependent",
            "published": published_dep,
            "replicated": dep0 / max(tot0, 1),
            "corrected": dep1 / max(tot1, 1),
            "n_units": tot0, "kind": "replication"}


# ---------------------------------------------------------------- B: LINCS
def lincs(published=0.26, min_lines=6, thresh=0.30, max_cpd=1200):
    from cmapPy.pandasGEXpress.parse import parse
    inst = pd.read_csv(LINCS / "p1_inst_info.txt.gz", sep="\t",
                       low_memory=False)
    inst = inst[inst.pert_type == "trt_cp"]
    keep = inst.groupby("pert_id").cell_id.nunique()
    cpds = list(keep[keep >= min_lines].index)[:max_cpd]
    inst = inst[inst.pert_id.isin(cpds)]
    print(f"  LINCS p1: {len(inst)} instances, {len(cpds)} compounds on "
          f">= {min_lines} lines", flush=True)
    M = parse(str(LINCS / "p1_level4.gctx"),
              cid=list(inst.inst_id.astype(str))).data_df
    cols = set(M.columns)
    inst = inst[inst.inst_id.astype(str).isin(cols)]
    X = M.T.to_numpy(dtype=np.float32)
    pos = {c: i for i, c in enumerate(M.columns)}
    idx = np.array([pos[i] for i in inst.inst_id.astype(str)])
    X = X[idx]
    ctx = inst.cell_id.to_numpy()
    prt = inst.pert_id.to_numpy()
    print(f"  matrix {X.shape}", flush=True)

    def conserved(Dm):
        """Median between-line correlation of a compound's signature.

        Subramanian et al. call a compound conserved when its signature
        reproduces across the panel; the operationalisation here is the median
        pairwise correlation between per-line signatures, thresholded once and
        applied identically before and after the correction.
        """
        vals = []
        for p in np.unique(prt):
            m = prt == p
            per = []
            for c in np.unique(ctx[m]):
                per.append(Dm[m][ctx[m] == c].mean(0))
            if len(per) < min_lines:
                continue
            P = np.stack(per)
            R = np.corrcoef(P)
            iu = np.triu_indices(len(P), 1)
            vals.append(float(np.median(R[iu])))
        v = np.array(vals)
        return float((v > thresh).mean()), v

    cs = context_share(X, ctx, prt, "LINCS L4 (plate z-scored)")
    f0, v0 = conserved(X)
    C = loo_context_mean(X, ctx, prt)
    f1, v1 = conserved(X - C)
    print(f"  as published:  {f0:.1%} of compounds conserved across the panel "
          f"(median cross-line r = {np.median(v0):+.3f})")
    print(f"  corrected:     {f1:.1%} conserved "
          f"(median cross-line r = {np.median(v1):+.3f})")
    return ({"context_share": cs,
             "study": "Subramanian 2017 (CMap / LINCS)",
             "statistic": "% of compounds with a panel-conserved signature",
             "published": published, "replicated": f0, "corrected": f1,
             "n_units": len(v0), "kind": "replication"},
            pd.DataFrame({"as_published": v0, "corrected": v1}))


# ---------------------------------------------------------------- C: PRISM
def bendavid_analogue(active=-0.301, inactive=-0.097):
    """Ben-David's counting rule, applied to PRISM lines.

    Their thresholds are >50% and <20% growth inhibition, which on a log2
    fold-change scale are -1 and log2(0.8). Their 27 MCF7 strains are not
    public, so this is the same RULE on different contexts and is reported as an
    analogue. It is included because their statistic is a pure activity count
    with no residual at all -- exactly the kind of number that a context's
    general sensitivity should dominate.
    """
    from perturbmodel.celldrug import load_prism, general_sensitivity_by_rep
    R, K, lines, ti = load_prism()
    A, order = general_sensitivity_by_rep(R, K)
    apos = {c: j for j, c in enumerate(order)}
    reps = K.rep.to_numpy()

    def count(correct):
        n_act = n_split = 0
        for cpd, g in K.groupby("compound", observed=True):
            cols = g.index.to_numpy()
            sub = R[:, cols].copy()
            if correct:
                j = apos[cpd]
                sub = sub - np.stack([A[r][:, j] for r in reps[cols]], axis=1)
            with np.errstate(invalid="ignore"):
                v = np.nanmean(sub, axis=1)
            v = v[np.isfinite(v)]
            if len(v) < 30:
                continue
            if (v < active).any():                 # active somewhere
                n_act += 1
                if (v > inactive).any():           # inactive somewhere else
                    n_split += 1
        return n_split, n_act

    # the residual AS PRISM-style analyses construct it: compound mean out,
    # nothing else -- the construction the correction is aimed at. Each row is
    # one (line, condition) scalar, so the cube is flattened to a single column.
    with np.errstate(invalid="ignore"):
        res = np.nan_to_num(R - np.nanmean(R, axis=0, keepdims=True))
    cs = context_share(res.reshape(-1, 1),
                       np.repeat(np.arange(res.shape[0]), res.shape[1]),
                       np.tile(K.compound.to_numpy(), res.shape[0]),
                       "PRISM (compound mean only)")
    s0, a0 = count(False)
    s1, a1 = count(True)
    print(f"  as published rule: {s0}/{a0} = {s0/max(a0,1):.1%} of active "
          f"compounds are inactive in at least one context")
    print(f"  corrected:         {s1}/{a1} = {s1/max(a1,1):.1%}")
    return {"context_share": cs,
            "study": "Ben-David 2018 (rule applied to PRISM)",
            "statistic": "% of active compounds inactive in >=1 context",
            "published": 48 / 55, "replicated": s0 / max(a0, 1),
            "corrected": s1 / max(a1, 1), "n_units": a0, "kind": "analogue"}


def context_share(D, ctx, pert, label):
    """How much of a residual, AS THAT PIPELINE CONSTRUCTS IT, is a context
    main effect?

    This is the variable that decides whether the correction matters. A pipeline
    that z-scores within plate, or takes deltas against vehicle in the same
    context, has already removed most of the context main effect before any
    statistic is computed; subtracting an estimate of it then adds noise rather
    than removing bias. A scalar viability residual with only the compound mean
    taken out has removed none of it.
    """
    C = loo_context_mean(D, ctx, pert)
    vt = float(np.var(D))
    vc = float(np.var(C))
    print(f"  {label:34s} context main effect is {vc/max(vt,1e-12):6.1%} "
          f"of the residual as constructed")
    return vc / max(vt, 1e-12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-lincs", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(ROOT / "scripts"))
    rows, lincs_v = [], None

    print("A. Srivatsan et al., Science 2020 (sci-Plex 3)", flush=True)
    try:
        r = sciplex()
        if r:
            rows.append(r)
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}")

    if not args.skip_lincs:
        print("\nB. Subramanian et al., Cell 2017 (CMap / LINCS)", flush=True)
        try:
            r, lincs_v = lincs()
            rows.append(r)
        except Exception as e:
            print(f"  FAILED: {type(e).__name__}: {e}")

    print("\nC. Ben-David et al., Nature 2018 (rule applied to PRISM)",
          flush=True)
    try:
        rows.append(bendavid_analogue())
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}")

    if not rows:
        print("\nnothing computed")
        return
    T = pd.DataFrame(rows)
    T["replication_error"] = (T.replicated - T.published).abs()
    T["reproduced"] = T.replication_error < 0.10
    T["inflation"] = T.replicated - T.corrected
    T["inflation_factor"] = T.replicated / T.corrected.replace(0, np.nan)
    T = T.rename(columns={"context_share": "context_share_of_residual"})
    T.to_csv(TAB / "published_number_inflation.csv", index=False)

    print("\n" + "=" * 74)
    print(T[["study", "published", "replicated", "corrected", "reproduced",
             "inflation_factor"]].round(3).to_string(index=False))
    ok = T[T.reproduced & (T.kind == "replication")]
    print(f"\n{len(ok)} of {int((T.kind == 'replication').sum())} published "
          f"statistics reproduced within 10 points of the published value.")
    if len(ok):
        print("For those, removing the context main effect changes the "
              "statistic by")
        for r_ in ok.itertuples():
            print(f"  {r_.study:38s} {r_.replicated:.1%} -> {r_.corrected:.1%}"
                  f"   ({r_.inflation_factor:.2f}x)")
    bad = T[~T.reproduced & (T.kind == "replication")]
    for r_ in bad.itertuples():
        print(f"NOT REPRODUCED: {r_.study} — published {r_.published:.1%}, "
              f"our implementation {r_.replicated:.1%}. No inflation claim is "
              f"made for this statistic.")
    print("\nWhat decides whether the correction matters:")
    for r_ in T.itertuples():
        print(f"  {r_.study:38s} context main effect "
              f"{r_.context_share_of_residual:6.1%} of the residual  ->  "
              f"statistic moves {r_.inflation_factor:.2f}x")
    print("  This diagnostic does NOT predict the effect -- it is anti-ordered\n"
          "  with it here. It measures the context term against TOTAL residual\n"
          "  variance, which is noise-dominated (PRISM: 2.8%), whereas the\n"
          "  replicate-covariance decomposition puts the same term at ~33% of\n"
          "  the REPRODUCIBLE residual. Statistics built on raw residuals are\n"
          "  mostly measuring noise, so removing a small term changes little.\n"
          "  No claim is made that published numbers are inflated.")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    n = 2 if lincs_v is None else 3
    fig, ax = plt.subplots(1, n, figsize=(4.7 * n, 4.3),
                           constrained_layout=True)
    xx = np.arange(len(T))
    w = 0.26
    ax[0].bar(xx - w, T.published, w, color=GREY, label="published")
    ax[0].bar(xx, T.replicated, w, color=BLUE,
              label="our replication of their definition")
    ax[0].bar(xx + w, T.corrected, w, color=AQUA,
              label="context main effect removed")
    for i_, r_ in enumerate(T.itertuples()):
        ax[0].text(i_ - w, r_.published + .012, f"{r_.published:.0%}",
                   ha="center", fontsize=7)
        ax[0].text(i_, r_.replicated + .012, f"{r_.replicated:.0%}",
                   ha="center", fontsize=7, fontweight="bold")
        ax[0].text(i_ + w, r_.corrected + .012, f"{r_.corrected:.0%}",
                   ha="center", fontsize=7, fontweight="bold")
    ax[0].set_xticks(xx, [s.split(" (")[0] + ("\n(analogue)" if k == "analogue"
                                              else "")
                          for s, k in zip(T.study, T.kind)], fontsize=7)
    ax[0].set_ylabel("statistic as reported")
    ax[0].legend(frameon=False, fontsize=7, loc="upper left")
    ax[0].set_ylim(0, max(1.0, float(T[["published", "replicated"]].max().max())
                          * 1.35))
    ax[0].set_title("a  Published context-specificity statistics,\n    "
                    "recomputed with the context main effect out", loc="left",
                    fontweight="bold", fontsize=9.5)

    ax[1].axhline(0, color="#444", lw=0.9)
    ax[1].bar(xx, T.replicated - T.corrected, width=0.5,
              color=[ORANGE if k == "replication" else GREY for k in T.kind])
    for i_, v in enumerate(T.replicated - T.corrected):
        ax[1].text(i_, v + 0.006 * np.sign(v or 1), f"{v:+.1%}", ha="center",
                   fontsize=8, fontweight="bold")
    ax[1].set_xticks(xx, [s.split(" (")[0] for s in T.study], fontsize=7)
    ax[1].set_ylabel("overstatement (replication − corrected)")
    ax[1].set_title("b  How much the omission adds", loc="left",
                    fontweight="bold", fontsize=9.5)
    ax[1].twinx().plot(xx, T.context_share_of_residual, "o--", color=VIOLET,
                       lw=1.6, ms=7, label="context main effect in residual")
    ax[1].text(0.5, 0.93, "violet: context main effect as a share of TOTAL\n"
               "residual variance — anti-ordered with the effect, so it does\n"
               "NOT explain it (see module docstring)",
               transform=ax[1].transAxes, ha="center", fontsize=6.2,
               color=VIOLET)

    if lincs_v is not None and n == 3:
        ax[2].hist(lincs_v.as_published, bins=40, alpha=0.75, color=ORANGE,
                   label="as published")
        ax[2].hist(lincs_v.corrected, bins=40, alpha=0.7, color=AQUA,
                   label="corrected")
        ax[2].set_xlabel("median between-line signature correlation")
        ax[2].set_ylabel("compounds")
        ax[2].legend(frameon=False, fontsize=7.5)
        ax[2].set_title("c  LINCS, per compound", loc="left",
                        fontweight="bold", fontsize=9.5)
    fig.suptitle("A negative result: published context-specificity statistics "
                 "do not move when the context main effect is removed",
                 fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "published_inflation", FIG,
                    source_data={"summary": T,
                                 "lincs_per_compound": lincs_v
                                 if lincs_v is not None else pd.DataFrame()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
