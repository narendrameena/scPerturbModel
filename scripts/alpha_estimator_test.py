#!/usr/bin/env python3
"""Does the line-effect correction remove a cell property, or eat the interaction?

Two analyses of the same Tahoe plates disagree sixty-fold about the context x
compound interaction: the scalar replicate analysis (RESULTS.md sec.31) reports
0.5% [0.0-1.5%] and the per-gene decomposition (sec.33) reports 34%. Two
candidate explanations were tested and failed -- dose pooling and gene set. This
script tests the third, which is the one that would matter most if true.

**The concern.** The correction subtracts each line's mean response across the
OTHER drugs. That is right when the thing being removed is general sensitivity: a
line that responds to everything. It is wrong when the interaction itself is
CORRELATED ACROSS DRUGS within a line -- a line unusually sensitive to kinase
inhibitors but not to DNA-damaging agents has a genuine context x compound
interaction, yet a mean over all ~95 drugs partially absorbs it. PRISM compound
residuals correlate at r-bar = +0.23, so this is not hypothetical. The correction
removes the rank-1 line component; whatever share of a real interaction is rank-1
in drug space goes with it.

Four estimators of the line effect, on identical data, identical pairs:

  none              no correction -- the pre-2026-09-03 convention, which
                    reports the cell property as interaction
  leave-one-drug    the current correction (sec.31)
  leave-one-class   exclude every drug sharing the query drug's annotated
                    mechanism, so a class-structured interaction cannot be
                    absorbed by its own class-mates
  half-split        estimate the line effect on a random half of drugs and
                    measure the interaction only on the OTHER half, so the
                    estimate and the measurement share no compound at all

If all three corrections agree, the correction is removing a cell property and
sec.31's null stands. If leave-one-class or half-split gives a substantially
larger interaction than leave-one-drug, the correction is over-aggressive and the
true value lies between -- which would make sec.31's 0.5% a lower bound rather
than an estimate.

Outputs: results/tables/alpha_estimator_test.csv
         figure bundle results/figures/00_manuscript/alpha_estimator_test/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    ap.add_argument("--n-boot", type=int, default=400)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    from perturbmodel.evaluation.delta_eval import (build_deltas,
                                                    load_pseudobulk,
                                                    responsive_genes)
    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond, keep_plate=True, exclude_plates=())
    resp = responsive_genes(DELTA, np.ones(len(G), bool))
    D = DELTA[:, resp].astype(np.float32)
    del DELTA
    K = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                      "drug": G.drug, "conc": G.conc, "plate": G.plate})
    print(f"{len(K)} conditions, {len(resp)} genes, "
          f"{K.drug.nunique()} drugs, {K.line.nunique()} lines", flush=True)

    md = pd.read_csv(TAB / "drug_metadata.csv")
    moa = dict(zip(md.drug.astype(str), md["moa-fine"].astype(str)))
    K["moa"] = [moa.get(d, "unclear") for d in K.drug]
    named = K[~K.moa.isin(["unclear", "nan", "", "None"])]
    print(f"{named.drug.nunique()} of {K.drug.nunique()} drugs carry a "
          f"mechanism annotation, {named.moa.nunique()} classes", flush=True)

    # shared response: leave-one-LINE-out mean per (drug, dose), exactly as in
    # the analysis under test
    resid, shared_num, shared_den = {}, 0.0, 0
    for (drug, conc), g in K.groupby(["drug", "conc"], observed=True):
        ii = g.i.to_numpy(); ln = g.line.to_numpy()
        if len(np.unique(ln)) < 2:
            continue
        tot = D[ii].sum(0)
        csum = {c: D[ii[ln == c]].sum(0) for c in np.unique(ln)}
        ccnt = {c: int((ln == c).sum()) for c in np.unique(ln)}
        for i, c in zip(ii, ln):
            n_out = len(ii) - ccnt[c]
            if n_out < 1:
                continue
            loo = (tot - csum[c]) / n_out
            resid[i] = D[i] - loo
            shared_num += float(np.mean(loo ** 2)); shared_den += 1
    shared = shared_num / shared_den
    print(f"shared component {shared:.5f} over {shared_den} conditions",
          flush=True)

    dv = K.drug.to_numpy(); mv = K.moa.to_numpy()

    def line_effect(exclude):
        """Per condition, the line's mean residual over drugs `exclude` allows.

        ``exclude(query_drug, other_drug) -> bool`` decides what is held out of
        each estimate. Estimated within (line, plate) because the delta is taken
        against controls on its own plate, then centred across lines within
        (drug, dose, plate) so what is subtracted is a CONTRAST between lines
        rather than the shared response wearing a line label.
        """
        alpha = {}
        for (ln, pl), g in K.groupby(["line", "plate"], observed=True):
            ii = [i for i in g.i.to_numpy() if i in resid]
            if len(ii) < 3:
                continue
            by = {}
            for i in ii:
                by.setdefault(dv[i], []).append(resid[i])
            by = {d: np.mean(v, axis=0) for d, v in by.items()}
            if len(by) < 3:
                continue
            for i in ii:
                keep = [v for d, v in by.items() if not exclude(i, d)]
                if len(keep) >= 2:
                    alpha[i] = np.mean(keep, axis=0)
        out = {}
        for (drug, conc, pl), g in K.groupby(["drug", "conc", "plate"],
                                             observed=True):
            ii = [i for i in g.i.to_numpy() if i in alpha]
            if len(ii) < 2:
                continue
            m = np.mean([alpha[i] for i in ii], axis=0)
            for i in ii:
                out[i] = alpha[i] - m
        return out

    moa_of = {i: mv[i] for i in range(len(K))}
    variants = {
        "none": None,
        "leave-one-drug": lambda i, d: d == dv[i],
        "leave-one-class": lambda i, d: (d == dv[i]) or
            (moa_of[i] not in ("unclear", "nan") and moa.get(d) == moa_of[i]),
    }

    def pairs_and_share(corr, restrict=None):
        """Cross-plate replicate covariance under one correction."""
        vals = []
        for (line, drug), g in K.groupby(["line", "drug"], observed=True):
            if restrict is not None and drug not in restrict:
                continue
            for conc, gd in g.groupby("conc", observed=True):
                ii = [i for i in gd.i.to_numpy() if i in resid]
                pl = [K.plate.iloc[i] for i in ii]
                for a in range(len(ii)):
                    for b in range(a + 1, len(ii)):
                        if pl[a] == pl[b]:
                            continue
                        ra, rb = resid[ii[a]], resid[ii[b]]
                        if corr is not None:
                            if ii[a] not in corr or ii[b] not in corr:
                                continue
                            ra = ra - corr[ii[a]]
                            rb = rb - corr[ii[b]]
                        vals.append(float(np.mean(ra * rb)))
        if not vals:
            return np.nan, np.nan, np.nan, 0
        v = np.array(vals)
        inter = float(v.mean())
        rng = np.random.default_rng(0)
        bs = [float(rng.choice(v, len(v), True).mean())
              for _ in range(args.n_boot)]
        def sh(x):
            return max(x, 0) / (shared + max(x, 0))
        return sh(inter), sh(float(np.percentile(bs, 2.5))), \
            sh(float(np.percentile(bs, 97.5))), len(v)

    rows = []
    for name, ex in variants.items():
        corr = line_effect(ex) if ex is not None else None
        s_, lo, hi, n = pairs_and_share(corr)
        rows.append({"estimator": name, "share": s_, "lo": lo, "hi": hi,
                     "n_pairs": n})
        print(f"  {name:18s} interaction share {s_:6.2%} "
              f"[{lo:.2%}-{hi:.2%}]   n = {n:,}", flush=True)

    # half-split: the line effect and the measurement share NO compound
    rng = np.random.default_rng(0)
    ud = np.array(sorted(K.drug.unique()))
    halves = []
    for rep_i in range(4):
        perm = rng.permutation(len(ud))
        hA = set(ud[perm[:len(ud) // 2]])
        corr = line_effect(lambda i, d: d in hA or d == dv[i])
        s_, lo, hi, n = pairs_and_share(corr, restrict=hA)
        halves.append(s_)
        print(f"  half-split rep {rep_i + 1}     interaction share {s_:6.2%}"
              f"   n = {n:,}", flush=True)
    rows.append({"estimator": "half-split", "share": float(np.mean(halves)),
                 "lo": float(np.min(halves)), "hi": float(np.max(halves)),
                 "n_pairs": n})

    T = pd.DataFrame(rows)
    T.to_csv(TAB / "alpha_estimator_test.csv", index=False)
    base = float(T[T.estimator == "leave-one-drug"].share.iloc[0])
    cls = float(T[T.estimator == "leave-one-class"].share.iloc[0])
    hs = float(T[T.estimator == "half-split"].share.iloc[0])
    print(f"\nVERDICT")
    print(f"  leave-one-drug (the current correction)  {base:.2%}")
    print(f"  leave-one-class                          {cls:.2%}")
    print(f"  half-split (no shared compound)          {hs:.2%}")
    if max(cls, hs) > 3 * max(base, 1e-4) and max(cls, hs) > 0.02:
        print("  The correction is OVER-AGGRESSIVE: holding out a drug's whole "
              "class, or\n  estimating the line effect on compounds that are "
              "not measured, recovers a\n  substantially larger interaction. "
              "sec.31's 0.5% is then a LOWER BOUND on\n  a class-structured "
              "interaction, not an estimate of it, and the headline\n  should be "
              "restated accordingly.")
    else:
        print("  The correction is NOT over-aggressive: excluding a drug's "
              "whole class, or\n  sharing no compound between the estimate and "
              "the measurement, leaves the\n  interaction essentially "
              "unchanged. What is removed behaves as a cell\n  property, and "
              "sec.31's null stands as an estimate.")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    xx = np.arange(len(T))
    cols = [GREY if e == "none" else ORANGE if e == "leave-one-drug" else VIOLET
            for e in T.estimator]
    ax[0].bar(xx, T.share, width=0.55, color=cols)
    ax[0].errorbar(xx, T.share, yerr=[T.share - T.lo, T.hi - T.share],
                   fmt="none", ecolor="#333", capsize=3, lw=1.1)
    for i_, r_ in enumerate(T.itertuples()):
        ax[0].text(i_, r_.hi + 0.004, f"{r_.share:.1%}", ha="center",
                   fontsize=9, fontweight="bold")
    ax[0].set_xticks(xx, [e.replace("-", "-\n", 1) for e in T.estimator],
                     fontsize=7)
    ax[0].set_ylabel("interaction share (Tahoe, same-dose replicates)")
    ax[0].set_title("a  How the line effect is estimated\n    decides the "
                    "answer", loc="left", fontweight="bold", fontsize=9.5)

    d = T[T.estimator != "none"]
    ax[1].bar(np.arange(len(d)), d.share, width=0.5,
              color=[ORANGE if e == "leave-one-drug" else VIOLET
                     for e in d.estimator])
    ax[1].set_xticks(np.arange(len(d)), [e for e in d.estimator], fontsize=7)
    ax[1].set_ylabel("interaction share")
    ax[1].text(0.5, 0.9, "if these agree, what is removed is a cell property;\n"
               "if the right-hand bars are larger, the correction is\n"
               "absorbing interaction that is correlated across drugs",
               transform=ax[1].transAxes, ha="center", fontsize=6.8,
               color="#444")
    ax[1].set_title("b  Is the correction eating the interaction?", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Testing whether the line-effect correction removes a cell "
                 "property or a class-structured interaction", fontsize=10.5,
                 x=0.005, ha="left", fontweight="bold")
    dd = save_figure(fig, "alpha_estimator_test", FIG,
                     source_data={"variants": T}, script=__file__)
    print(f"figure bundle -> {dd}")


if __name__ == "__main__":
    main()
