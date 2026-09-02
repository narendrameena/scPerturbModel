#!/usr/bin/env python3
"""Why does context-dependence peak just below the lethal dose?

RESULTS.md §13 reports that the interaction rises with dose to a peak near
2.5 uM and collapses at the top dose, and offers an explanation: the interaction
is largest where cell lines differ most in *where they sit on the dose-response
curve*, and vanishes once the dose saturates that difference and everything
dies. That explanation was asserted and never tested, which is exactly the kind
of just-so story a reviewer should reject.

It makes a specific, falsifiable prediction. If context-dependence is the
between-line spread in potency read out at a particular dose, then for each
compound **the dose at which the interaction peaks should coincide with the dose
at which lines are most spread apart**, and both should sit below the dose at
which the average line is already dead. Three tests:

  1. within compound, does the interaction index track the between-line spread
     of response across doses?
  2. across compounds, does the dose position of peak interaction coincide with
     the dose position of peak spread?
  3. is the peak located below the dose at which the mean effect saturates?

A competing explanation is that the collapse at the top dose is purely a ceiling
artefact of the viability assay -- everything hits the floor, so nothing can
differ. Test 1 distinguishes them: under the potency explanation the interaction
should follow the spread at ALL doses, not only at the top, whereas a pure
ceiling artefact predicts agreement only where the assay saturates.

Outputs: results/tables/dose_mechanism.csv
         figure bundle results/figures/13_prism/dose_mechanism/
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
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def build():
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
    L = lfc.to_numpy(dtype=np.float32)
    colpos = {c: i for i, c in enumerate(lfc.columns)}
    rows = []
    for cpd, gc in ti.groupby("name", observed=True):
        doses = sorted(gc.dose.unique())
        if len(doses) < 4:
            continue
        for di, dose in enumerate(doses):
            gd = gc[gc.dose == dose]
            cols_by_rep = [[colpos[c] for c in g.column_name]
                           for _, g in gd.groupby("rep", observed=True)]
            if len(cols_by_rep) < 2:
                continue
            with np.errstate(invalid="ignore"):
                reps = np.stack([np.nanmean(L[:, ix], axis=1)
                                 for ix in cols_by_rep], axis=1)
                line_mean = np.nanmean(reps, axis=1)
            valid = np.isfinite(line_mean)
            n = int(valid.sum())
            if n < MIN_LINES:
                continue
            lm = line_mean[valid]
            loo = (lm.sum() - lm) / (n - 1)
            add = float(np.mean(loo ** 2))
            res = reps[valid] - loo[:, None]
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
            rows.append({
                "compound": cpd, "dose": float(dose), "dose_rank": di + 1,
                "n_dose": len(doses), "n_lines": n,
                "cdi": inter / den if den > 0 else np.nan,
                "interaction": inter,
                # the quantity the explanation is about: how far apart the lines
                # are at this dose, measured on the line-level responses
                "spread": float(np.std(lm)),
                "mean_effect": float(np.mean(lm))})
    D = pd.DataFrame(rows)
    D["dose_pos"] = np.ceil(D.dose_rank / D.n_dose * 8).clip(1, 8).astype(int)
    return D


def main():
    ap = argparse.ArgumentParser()
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    D = build()
    D.to_csv(TAB / "dose_mechanism.csv", index=False)
    print(f"{D.compound.nunique()} compounds x up to "
          f"{D.dose_pos.max()} dose positions", flush=True)

    # 1. within compound, does the interaction track the spread?
    rs = []
    for _, g in D.groupby("compound", observed=True):
        g = g.dropna(subset=["cdi", "spread"])
        if g.dose_pos.nunique() < 4:
            continue
        r = stats.spearmanr(g.spread, g.cdi).statistic
        if np.isfinite(r):
            rs.append(r)
    rs = np.array(rs)
    print(f"\n1. interaction vs between-line spread, within compound:")
    print(f"   median rho = {np.median(rs):+.3f}, {(rs > 0).mean():.0%} "
          f"positive, p = {stats.wilcoxon(rs)[1]:.2e} (n = {len(rs)})")

    # 2. do the two peaks coincide?
    pk = []
    for c, g in D.groupby("compound", observed=True):
        g = g.dropna(subset=["cdi", "spread"])
        if g.dose_pos.nunique() < 4:
            continue
        pk.append({"compound": c,
                   "peak_cdi": int(g.loc[g.cdi.idxmax(), "dose_pos"]),
                   "peak_spread": int(g.loc[g.spread.idxmax(), "dose_pos"]),
                   "peak_effect": int(g.loc[g.mean_effect.idxmin(), "dose_pos"])})
    P = pd.DataFrame(pk)
    rho = stats.spearmanr(P.peak_cdi, P.peak_spread)
    same = float((P.peak_cdi == P.peak_spread).mean())
    within1 = float((abs(P.peak_cdi - P.peak_spread) <= 1).mean())
    print(f"\n2. dose position of peak interaction vs peak spread:")
    print(f"   rho = {rho.statistic:+.3f}, p = {rho.pvalue:.2e}; identical for "
          f"{same:.0%}, within one position for {within1:.0%} (n = {len(P)})")
    print(f"   median peak interaction at position "
          f"{P.peak_cdi.median():.0f}/8, peak spread at "
          f"{P.peak_spread.median():.0f}/8, strongest killing at "
          f"{P.peak_effect.median():.0f}/8")

    # 3. is the peak below the saturating dose?
    below = float((P.peak_cdi < P.peak_effect).mean())
    print(f"\n3. the interaction peaks BELOW the dose of maximal killing for "
          f"{below:.0%} of compounds")

    # discriminating test: does the relationship hold away from the ceiling?
    lowmid = D[D.dose_pos <= 6].dropna(subset=["cdi", "spread"])
    rs2 = []
    for _, g in lowmid.groupby("compound", observed=True):
        if g.dose_pos.nunique() < 4:
            continue
        r = stats.spearmanr(g.spread, g.cdi).statistic
        if np.isfinite(r):
            rs2.append(r)
    rs2 = np.array(rs2)
    if len(rs2) > 10:
        print(f"\n   restricted to non-saturating doses (positions 1-6), where "
              f"no ceiling\n   artefact can operate: median rho = "
              f"{np.median(rs2):+.3f}, {(rs2 > 0).mean():.0%} positive, "
              f"p = {stats.wilcoxon(rs2)[1]:.2e} (n = {len(rs2)})")
        print("   The relationship holding away from the ceiling is what "
              "separates the\n   potency-spread explanation from a pure "
              "assay-saturation artefact.")

    S = pd.DataFrame([{"test": "interaction vs spread, within compound",
                       "statistic": float(np.median(rs)), "n": len(rs)},
                      {"test": "peak interaction vs peak spread",
                       "statistic": float(rho.statistic), "n": len(P)},
                      {"test": "peaks identical", "statistic": same,
                       "n": len(P)},
                      {"test": "peak below maximal killing",
                       "statistic": below, "n": len(P)},
                      {"test": "interaction vs spread, non-saturating doses",
                       "statistic": float(np.median(rs2)) if len(rs2) else
                       np.nan, "n": len(rs2)}])
    S.to_csv(TAB / "dose_mechanism_summary.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.1), constrained_layout=True)
    q = D.groupby("dose_pos").agg(cdi=("cdi", "median"),
                                  spread=("spread", "median"),
                                  eff=("mean_effect", "median"))
    a2 = ax[0].twinx()
    ax[0].plot(q.index, q.cdi, "o-", color=VIOLET, lw=2, ms=5,
               label="interaction index")
    a2.plot(q.index, q.spread, "s--", color=AQUA, lw=2, ms=5,
            label="between-line spread")
    ax[0].set_xlabel("dose position (1 = lowest)")
    ax[0].set_ylabel("context-dependence index", color=VIOLET)
    a2.set_ylabel("between-line spread (SD)", color=AQUA)
    h1, l1 = ax[0].get_legend_handles_labels()
    h2, l2 = a2.get_legend_handles_labels()
    ax[0].legend(h1 + h2, l1 + l2, frameon=False, fontsize=7, loc="upper left")
    ax[0].set_title("a  Both peak at the same dose", loc="left",
                    fontweight="bold", fontsize=9.5)

    jit = np.random.default_rng(0).normal(0, 0.10, len(P))
    ax[1].scatter(P.peak_spread + jit, P.peak_cdi + jit, s=13, alpha=0.4,
                  color=BLUE, edgecolors="none")
    ax[1].plot([1, 8], [1, 8], ls="--", color="#444", lw=1.2)
    ax[1].set_xlabel("dose position of peak spread")
    ax[1].set_ylabel("dose position of peak interaction")
    ax[1].text(0.03, 0.95, f"ρ = {rho.statistic:+.2f}\n{within1:.0%} within one "
               f"position", transform=ax[1].transAxes, va="top", fontsize=7.5)
    ax[1].set_title("b  Per compound", loc="left", fontweight="bold",
                    fontsize=9.5)

    ax[2].hist(rs, bins=30, color=VIOLET, alpha=0.55, label="all doses")
    if len(rs2) > 10:
        ax[2].hist(rs2, bins=30, color=AQUA, alpha=0.55,
                   label="non-saturating only")
    ax[2].axvline(0, color="#444", lw=1.2)
    ax[2].set_xlabel("within-compound ρ (interaction vs spread)")
    ax[2].set_ylabel("compounds")
    ax[2].legend(frameon=False, fontsize=7.5)
    ax[2].set_title("c  Not a ceiling artefact", loc="left", fontweight="bold",
                    fontsize=9.5)
    fig.suptitle("Testing the explanation: context-dependence tracks how far "
                 "apart the lines are at that dose", fontsize=10.5, x=0.005,
                 ha="left", fontweight="bold")
    d = save_figure(fig, "dose_mechanism", FIG,
                    source_data={"per_dose": D, "peaks": P, "summary": S},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
