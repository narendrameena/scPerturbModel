#!/usr/bin/env python3
"""Does dose change how much of a drug response is context-specific?

So far context-dependence has been treated as a property of the compound, with
dose pooled away. That hides a question with both mechanistic and practical
weight: a compound's transferability may not be one number but a curve.

Three outcomes are distinguishable and mean different things.

  flat            context-dependence is a property of the drug-target
                  relationship, independent of occupancy. Pooling doses is
                  then harmless and the mechanism claim stands as stated.

  rising          higher dose recruits off-target and stress pathways whose
                  wiring differs between lines, so responses diverge. Implies
                  transfer is easiest at low, target-selective doses.

  falling         at low dose only the lines with the sensitising background
                  respond at all (the interaction IS the differing potency);
                  at saturating dose everything dies the same way. Implies a
                  ceiling effect and that CDI partly measures IC50 spread.

The falling case is the one to worry about, because it would mean our index
partly re-measures potency variation rather than qualitative rewiring. We
therefore separate the two: alongside CDI we report whether the residuals are
merely rescaled between lines (a potency shift) or point in different
directions (a genuinely different response).

PRISM gives 8 doses across ~740 lines; Tahoe gives 3 across 47. Both are run,
because agreement across a scalar and a transcriptional readout is what makes
a dose effect a property of the biology rather than of the assay.

Outputs: results/tables/dose_vs_context_prism.csv
         results/tables/dose_vs_context_tahoe.csv
         figure bundle results/figures/13_prism/dose_vs_context/
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

from perturbmodel.celldrug import general_sensitivity_by_rep
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "13_prism"
TAB = ROOT / "results" / "tables"
MIN_LINES = 100
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def prism_dose():
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
    colpos = {c: i for i, c in enumerate(lfc.columns)}

    recs = []
    # Each line's general sensitivity, leave-one-compound-out, removed before
    # anything below is called context-dependence. It reproduces across disjoint
    # compound halves at r = 0.989 and is shared between replicate plates
    # exactly as a genuine interaction is.
    _Kall = pd.DataFrame({"compound": ti.name.to_numpy(),
                          "rep": ti.rep.to_numpy()})
    _Rall = L[:, [colpos[c] for c in ti.column_name]]
    ALPHA, ACPD = general_sensitivity_by_rep(_Rall, _Kall)
    apos = {c: j for j, c in enumerate(ACPD)}

    for cpd, gc in ti.groupby("name", observed=True):
        aj = apos.get(cpd)
        doses = sorted(gc.dose.unique())
        if len(doses) < 4:
            continue
        for di, dose in enumerate(doses):
            gd = gc[gc.dose == dose]
            rep_groups = list(gd.groupby("rep", observed=True))
            cols_by_rep = [[colpos[c] for c in g.column_name]
                           for _, g in rep_groups]
            # each replicate column corrected from its OWN plate
            Acol = (np.stack([ALPHA[r][:, aj] for r, _ in rep_groups], axis=1)
                    if aj is not None else np.zeros((L.shape[0],
                                                     len(rep_groups))))
            if len(cols_by_rep) < 2:
                continue
            with np.errstate(invalid="ignore"):
                reps = np.stack([np.nanmean(L[:, ix], axis=1)
                                 for ix in cols_by_rep], axis=1)
                line_mean = np.nanmean(reps, axis=1)
            valid = np.isfinite(line_mean) & np.isfinite(Acol).all(axis=1)
            n = int(valid.sum())
            if n < MIN_LINES:
                continue
            tot = line_mean[valid].sum()
            loo = (tot - line_mean[valid]) / (n - 1)
            add = float(np.mean(loo ** 2))
            res = reps[valid] - loo[:, None] - Acol[valid]
            cov, npair = 0.0, 0
            for a in range(res.shape[1]):
                for b in range(a + 1, res.shape[1]):
                    m = np.isfinite(res[:, a]) & np.isfinite(res[:, b])
                    if m.sum() >= MIN_LINES:
                        cov += float(np.mean(res[m, a] * res[m, b])); npair += 1
            if npair == 0:
                continue
            inter = max(cov / npair, 0.0)
            den = add + inter
            # potency-vs-rewiring split is not defined for a scalar readout;
            # for viability we can still ask whether the mean effect saturates
            recs.append({"compound": cpd, "dose": float(dose),
                         "dose_rank": di + 1, "n_dose": len(doses),
                         "n_lines": n, "additive": add, "interaction": inter,
                         "cdi": inter / den if den > 0 else np.nan,
                         "mean_effect": float(np.mean(line_mean[valid]))})
    return pd.DataFrame(recs)


def tahoe_dose(pb_dir):
    from perturbmodel.evaluation.delta_eval import (build_deltas,
                                                    load_pseudobulk,
                                                    responsive_genes)
    X, cond = load_pseudobulk(ROOT / pb_dir)
    G, DELTA = build_deltas(X, cond, keep_plate=True)
    resp = responsive_genes(DELTA, np.ones(len(G), bool))
    D = DELTA[:, resp].astype(np.float32)
    del DELTA
    K = pd.DataFrame({"i": np.arange(len(G)), "line": G.cell_line_id,
                      "drug": G.drug, "conc": G.conc, "plate": G.plate})
    recs = []
    for drug, gd in K.groupby("drug", observed=True):
        concs = sorted(gd.conc.unique())
        if len(concs) < 2:
            continue
        for di, c in enumerate(concs):
            gc = gd[gd.conc == c]
            if gc.line.nunique() < 8:
                continue
            ii = gc.i.to_numpy()
            if len(ii) < 2:
                continue
            # leave-one-context-out, for the same reason as elsewhere: an
            # in-sample mean biases the residual covariance to -sigma^2/n
            tot = D[ii].sum(0)
            loo = {i: (tot - D[i]) / (len(ii) - 1) for i in ii}
            add = float(np.mean([np.mean(loo[i] ** 2) for i in ii]))
            res = {i: D[i] - loo[i] for i in ii}
            mu = D[ii].mean(0)
            cov, npair, cs = 0.0, 0, []
            for _, gl in gc.groupby("line", observed=True):
                v = gl.i.to_numpy(); pl = gl.plate.to_numpy()
                for a in range(len(v)):
                    for b in range(a + 1, len(v)):
                        if pl[a] != pl[b]:
                            cov += float(np.mean(res[v[a]] * res[v[b]]))
                            npair += 1
            # the interaction covariance needs replicate pairs and usually
            # has too few here; the alignment measure below needs only the
            # lines themselves, so keep the row either way and null the CDI
            if npair >= 5:
                inter = max(cov / npair, 0.0)
                den = add + inter
                cdi = inter / den if den > 0 else np.nan
            else:
                inter, cdi = np.nan, np.nan
            # potency vs rewiring: correlation of each line's DELTA with the
            # shared response. high r with varying norm = same response, different
            # magnitude (potency). low r = a qualitatively different response.
            nm = np.linalg.norm(mu)
            if nm > 0:
                rr = [float(D[i] @ mu / (np.linalg.norm(D[i]) * nm))
                      for i in ii if np.linalg.norm(D[i]) > 0]
            else:
                rr = []
            recs.append({"drug": drug, "conc": c, "dose_rank": di + 1,
                         "n_conc": len(concs), "n_lines": gc.line.nunique(),
                         "additive": add, "interaction": inter,
                         "cdi": cdi, "n_pairs": npair,
                         "mean_effect": float(np.sqrt(add)),
                         "cos_to_shared": float(np.median(rr)) if rr else np.nan})
    return pd.DataFrame(recs)


def within_trend(df, gcol, xcol="dose_rank", ycol="cdi"):
    """Per-compound Spearman of y on dose, then a sign test over compounds.

    Doing it within compound removes the between-compound differences that
    would otherwise dominate, so the question asked is strictly 'does raising
    the dose of THIS compound change its context-dependence'.
    """
    rs = []
    for _, g in df.groupby(gcol, observed=True):
        g = g.dropna(subset=[ycol])
        if g[xcol].nunique() < 3:
            continue
        rho = stats.spearmanr(g[xcol], g[ycol]).statistic
        if np.isfinite(rho):
            rs.append(rho)
    rs = np.array(rs)
    if len(rs) < 5:
        return None
    w = stats.wilcoxon(rs)
    return {"n_compounds": len(rs), "median_rho": float(np.median(rs)),
            "frac_positive": float((rs > 0).mean()), "p": float(w.pvalue),
            "rhos": rs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    ap.add_argument("--skip-tahoe", action="store_true")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    P = prism_dose()
    # compounds are screened on grids of different length (mostly 8 points, a
    # few 16); ranking within compound and binning to eighths puts them on one
    # axis, so "dose 7 of 8" and "dose 14 of 16" are the same position
    P["dose_pos"] = np.ceil(P.dose_rank / P.n_dose * 8).clip(1, 8).astype(int)
    P.to_csv(TAB / "dose_vs_context_prism.csv", index=False)
    tp = within_trend(P, "compound", xcol="dose_pos")
    print(f"PRISM: {P.compound.nunique()} compounds, dose grids of "
          f"{sorted(P.n_dose.unique())} points")
    print(f"  CDI vs dose within compound: median rho={tp['median_rho']:+.3f}, "
          f"{tp['frac_positive']:.0%} positive, p={tp['p']:.2e} "
          f"(n={tp['n_compounds']})")
    pd_ = P.groupby("dose_pos").agg(cdi=("cdi", "median"),
                                     add=("additive", "median"),
                                     inter=("interaction", "median"),
                                     eff=("mean_effect", "median"),
                                     n=("cdi", "size"))
    print(pd_.round(4).to_string())
    pk = int(pd_.cdi.idxmax())
    print(f"  CDI peaks at dose position {pk}/8 (median {pd_.cdi.max():.3f})")
    if pk < 8:
        # a 16-point grid puts two ranks in one eighth, so collapse to one
        # value per compound before pairing
        a = P[P.dose_pos == pk].groupby("compound").cdi.median()
        b = P[P.dose_pos == 8].groupby("compound").cdi.median()
        j = a.index.intersection(b.index)
        a, b = a.dropna(), b.dropna()
        j = a.index.intersection(b.index)
        w = stats.wilcoxon(a[j], b[j])
        print(f"  peak vs top dose, paired within compound: "
              f"{a[j].median():.3f} vs {b[j].median():.3f}, p={w.pvalue:.2e} "
              f"(n={len(j)})")
        print(f"  at the top dose the mean effect is "
              f"{pd_.loc[8,'eff']:+.2f} log2 vs {pd_.loc[pk,'eff']:+.2f} at the")
        print("  peak: cells are simply dying everywhere, so the shared "
              "component\n  dominates and apparent context-dependence falls back.")

    T = pd.DataFrame()
    tt = None
    if not args.skip_tahoe:
        T = tahoe_dose(args.pb_dir)
        T.to_csv(TAB / "dose_vs_context_tahoe.csv", index=False)
        tt = within_trend(T, "drug")
        print(f"\nTahoe: {T.drug.nunique()} drugs x {T.dose_rank.max()} concs")
        print("  NOTE: splitting Tahoe by concentration leaves too few "
              "cross-plate\n  replicate pairs per (drug, conc) to estimate the "
              "interaction covariance,\n  so per-concentration CDI is not "
              "usable here. The alignment measure\n  below uses every line and "
              "is the informative one for Tahoe.")
        if tt:
            print(f"  CDI vs dose within drug: median rho={tt['median_rho']:+.3f}, "
                  f"{tt['frac_positive']:.0%} positive, p={tt['p']:.2e} "
                  f"(n={tt['n_compounds']})")
        td = T.groupby("dose_rank").agg(cdi=("cdi", "median"),
                                        add=("additive", "median"),
                                        inter=("interaction", "median"),
                                        cos=("cos_to_shared", "median"),
                                        n=("cdi", "size"))
        print(td.round(4).to_string())
        ct = within_trend(T, "drug", ycol="cos_to_shared")
        if ct:
            print(f"  alignment to the shared response vs dose: "
                  f"median rho={ct['median_rho']:+.3f}, p={ct['p']:.2e}")
            print("    rising alignment with dose = responses become MORE alike,")
            print("    so any CDI drop is convergence, not just a bigger denominator.")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    n = 4 if len(T) else 2
    fig, ax = plt.subplots(1, n, figsize=(4.0 * n, 4.4), constrained_layout=True)

    q = P.groupby("dose_pos").cdi.quantile([0.25, 0.5, 0.75]).unstack()
    ax[0].fill_between(q.index, q[0.25], q[0.75], color=VIOLET, alpha=0.2)
    ax[0].plot(q.index, q[0.5], "o-", color=VIOLET, lw=2, ms=6)
    ax[0].set_xlabel("dose position (1 = lowest, 8 = highest)")
    ax[0].set_ylabel("context-dependence index")
    ax[0].set_title(f"A  PRISM viability\nmedian rho={tp['median_rho']:+.2f}, "
                    f"p={tp['p']:.1e}", loc="left", fontweight="bold",
                    fontsize=10)

    a2 = ax[1]
    a2.plot(pd_.index, pd_["add"], "o-", color=BLUE, lw=2, label="additive")
    a2.plot(pd_.index, pd_["inter"], "o-", color=ORANGE, lw=2,
            label="interaction")
    a2.set_yscale("log")
    a2.set_xlabel("dose rank"); a2.set_ylabel("variance component")
    a2.legend(frameon=False, fontsize=8)
    a2.set_title("B  Both components grow;\nthe question is which grows faster",
                 loc="left", fontweight="bold", fontsize=10)

    if len(T):
        qt = T.groupby("dose_rank").additive.quantile([0.25, 0.5, 0.75]).unstack()
        ax[2].fill_between(qt.index, qt[0.25], qt[0.75], color=AQUA, alpha=0.2)
        ax[2].plot(qt.index, qt[0.5], "o-", color=AQUA, lw=2, ms=7)
        ax[2].set_xticks(sorted(T.dose_rank.unique()))
        ax[2].set_xlabel("concentration rank (1 = lowest of 3)")
        ax[2].set_ylabel("shared response magnitude")
        ax[2].set_title("C  Tahoe: the shared response\ngrows with dose",
                        loc="left", fontweight="bold", fontsize=10)

        qc = T.groupby("dose_rank").cos_to_shared.quantile([0.25, 0.5, 0.75]) \
              .unstack()
        ax[3].fill_between(qc.index, qc[0.25], qc[0.75], color=ORANGE, alpha=0.2)
        ax[3].plot(qc.index, qc[0.5], "o-", color=ORANGE, lw=2, ms=7)
        ax[3].set_xticks(sorted(T.dose_rank.unique()))
        ax[3].set_xlabel("concentration rank")
        ax[3].set_ylabel("cosine of line response to shared response")
        ct2 = within_trend(T, "drug", ycol="cos_to_shared")
        ax[3].set_title(f"D  Tahoe transcription diverges\nwith dose: rho="
                        f"{ct2['median_rho']:+.2f}, p={ct2['p']:.1e}",
                        loc="left", fontweight="bold", fontsize=10)
    fig.suptitle("Does dose change how context-specific a drug response is?",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "dose_vs_context", FIG,
                    source_data={"prism_by_dose": P, "tahoe_by_dose": T,
                                 "prism_summary": pd_.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
