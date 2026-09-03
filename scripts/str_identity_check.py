#!/usr/bin/env python3
"""STR profiling: an orthogonal identity check we had overlooked.

§18 concluded that cross-laboratory transfer improves when cell-line identity is
verified from drug-response fingerprints, and §23 walked that back in part,
noting that fingerprinting cannot separate genuine culture divergence from
metric noise and that no molecular identity check was available across these
datasets.

That last statement was wrong. PRISM ships `passed_str_profiling` for every cell
line, and **248 of 1,120 lines (22%) fail it**. Short-tandem-repeat profiling is
the field's standard molecular test of cell-line identity -- it is what
authentication services run, and failure means the culture does not match the
reference genotype for the line it is labelled as. It is entirely independent of
drug response, which makes it the orthogonal check the identity argument needed.

Three questions:

  1. do STR-failing lines transfer worse between laboratories than STR-passing
     ones? If identity is what limits transfer, they should.
  2. are STR failures enriched among the lines whose response fingerprint fails
     to identify them? Two independent identity signals agreeing would be strong
     evidence; disagreement would say the fingerprint measures something else.
  3. does excluding STR failures recover cross-laboratory agreement, and how
     does that compare with fingerprint-based selection?

Question 1 is the one that matters, because STR status is assigned without any
reference to the response data and cannot be circular.

Outputs: results/tables/str_identity_check.csv
         figure bundle results/figures/16_crosslab/str_identity/
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
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "16_crosslab"
TAB = ROOT / "results" / "tables"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    ap = argparse.ArgumentParser()
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(ROOT / "scripts"))
    import cross_lab_reproducibility as clr

    ci = pd.read_csv(PR / "secondary-screen-cell-line-info.csv")
    # The STR status cannot be joined on depmap_id: every one of the 248
    # STR-failing rows has depmap_id = NaN, so that join returned 747 True and
    # 1 False and the whole test was vacuous -- which is how the published claim
    # "all 738 lines with response data pass STR profiling" came about. The
    # status is only recoverable from row_name (PR500_ACH-000010_FAILED_STR).
    rn = ci.row_name.astype(str)
    ach = rn.str.extract(r"(ACH-\d+)")[0]
    failed = rn.str.contains("FAILED_STR", na=False)
    str_ok = {}
    for a, f in zip(ach, failed):
        if isinstance(a, str):
            str_ok[a] = (not f) and str_ok.get(a, True)
    n_fail_lines = sum(1 for v in str_ok.values() if not v)
    print(f"STR status parsed from row_name: {len(str_ok)} lines, "
          f"{n_fail_lines} with at least one FAILED_STR record", flush=True)
    print(f"  (the depmap_id column carries the status for none of the "
          f"failures)", flush=True)

    print("building residuals ...", flush=True)
    p_full, _, _ = clr.prism_residuals()
    G = clr.gdsc_residuals()
    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(r"DepMap;\s*(ACH-\d+)")[0]
    cos = {}
    for d, xr in zip(dep, smp.cross_references.astype(str)):
        if isinstance(d, str):
            for cid in re.findall(r"Cosmic(?:-CLP)?;\s*(\d+)", xr):
                cos.setdefault(int(cid), d)
    gb = {}
    for tag in ("GDSC1", "GDSC2"):
        g = G[tag].copy()
        g["dep"] = g.COSMIC_ID.map(cos)
        g = g[g.dep.notna()]
        for k, gg in g.groupby("k", observed=True):
            gb.setdefault(k, []).append(gg.groupby("dep").Z_SCORE.mean())
    gboth = {k: pd.concat(v).groupby(level=0).mean() for k, v in gb.items()}
    P = {norm(k): v for k, v in p_full.items()}

    # ---- Q1: does STR status predict cross-laboratory agreement? ----------
    # Correlate PRISM and GDSC residuals ACROSS COMPOUNDS for each line, so the
    # unit of analysis is the cell line and STR status can be tested on it.
    shared = sorted(set(P) & set(gboth))
    per_line = {}
    for cpd in shared:
        a, b = P[cpd], gboth[cpd]
        for ln in a.index.intersection(b.index):
            per_line.setdefault(ln, []).append((float(a[ln]), float(b[ln])))
    rows = []
    for ln, pairs in per_line.items():
        if len(pairs) < 25:
            continue
        u = np.array([p[0] for p in pairs]); v = np.array([p[1] for p in pairs])
        if np.std(u) == 0 or np.std(v) == 0:
            continue
        rows.append({"line": ln, "n_compounds": len(pairs),
                     "cross_lab_r": float(stats.spearmanr(u, v).statistic),
                     "str_pass": str_ok.get(ln, np.nan)})
    L = pd.DataFrame(rows)
    L["str_pass"] = L.str_pass.map({True: True, False: False})
    ok = L[L.str_pass == True].cross_lab_r
    bad = L[L.str_pass == False].cross_lab_r
    print(f"\n1. cross-laboratory agreement per cell line "
          f"({len(L)} lines, >=25 shared compounds each):")
    print(f"   STR pass  n={len(ok):4d}  median r = {ok.median():+.3f}")
    print(f"   STR fail  n={len(bad):4d}  median r = {bad.median():+.3f}")
    if len(ok) > 5 and len(bad) > 5:
        u_, p_ = stats.mannwhitneyu(ok, bad, alternative="two-sided")
        print(f"   difference {ok.median()-bad.median():+.3f}, "
              f"Mann-Whitney p = {p_:.3e}")
        print("   STR status is assigned without reference to drug response, so "
              "this cannot\n   be circular the way fingerprint selection can.")
    else:
        p_ = np.nan

    # ---- Q2: do the two identity signals agree? ---------------------------
    idf = TAB / "cross_lab_identity_viability.csv"
    agree = None
    if idf.exists():
        ID = pd.read_csv(idf)
        ID["str_pass"] = ID.line.map(str_ok)
        m = ID[ID.str_pass.isin([True, False])]
        if len(m) > 20:
            tab = pd.crosstab(m.str_pass, m.reciprocal_best)
            print(f"\n2. agreement between the two identity signals:")
            print(tab.to_string())
            if tab.shape == (2, 2):
                odds, pf = stats.fisher_exact(tab.to_numpy())
                print(f"   Fisher odds ratio {odds:.2f}, p = {pf:.3f}")
                print("   An odds ratio near 1 means the response fingerprint "
                      "and STR profiling are\n   measuring different things — "
                      "which is consistent with §23's conclusion that\n   "
                      "best-hit failure is largely metric noise.")
            agree = tab
            rk = m.dropna(subset=["rank_of_id_match"])
            if len(rk) > 20:
                r1 = rk[rk.str_pass == True].rank_of_id_match
                r0 = rk[rk.str_pass == False].rank_of_id_match
                if len(r0) > 5:
                    print(f"   median fingerprint rank: STR pass "
                          f"{r1.median():.0f}, STR fail {r0.median():.0f} "
                          f"(p = {stats.mannwhitneyu(r1, r0)[1]:.3f})")

    # ---- Q3: how much does excluding STR failures recover? ----------------
    def restrict(store, keep):
        return {k: v[v.index.isin(keep)] for k, v in store.items()}

    ceil = np.sqrt(0.473 * 0.438)
    keep_all = set(L.line)
    keep_str = set(L[L.str_pass == True].line)
    res = []
    for lab, keep in (("all lines", keep_all), ("STR-passing only", keep_str)):
        rr = clr.per_compound_corr(restrict(P, keep), restrict(gboth, keep), lab)
        med = float(np.median([x["rho"] for x in rr])) if rr else np.nan
        res.append({"selection": lab, "n_lines": len(keep),
                    "n_compounds": len(rr), "median_rho": med,
                    "frac_of_ceiling": med / ceil})
        print(f"\n3. {lab:20s} {len(keep):4d} lines  r={med:.3f}  "
              f"{med/ceil:.0%} of ceiling")
    S = pd.DataFrame(res)
    out = pd.concat([L.assign(analysis="per_line"),
                     S.assign(analysis="selection")], ignore_index=True)
    out.to_csv(TAB / "str_identity_check.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.1), constrained_layout=True)
    d = [ok.dropna(), bad.dropna()]
    bp = ax[0].boxplot(d, showfliers=False, patch_artist=True,
                       medianprops=dict(color="black", lw=1.5))
    for pch, c in zip(bp["boxes"], [AQUA, ORANGE]):
        pch.set_facecolor(c); pch.set_alpha(0.8)
    ax[0].set_xticks([1, 2], [f"STR pass\n(n={len(ok)})",
                              f"STR fail\n(n={len(bad)})"], fontsize=8)
    for i_, dd in enumerate(d):
        if len(dd):
            ax[0].text(i_ + 1, np.median(dd) + 0.02, f"{np.median(dd):.3f}",
                       ha="center", fontsize=8.5, fontweight="bold")
    ax[0].axhline(0, color="#444", lw=0.9)
    ax[0].set_ylabel("cross-laboratory r for that cell line")
    ax[0].set_title(f"a  STR status vs transfer\np = {p_:.1e}" if
                    np.isfinite(p_) else "a  STR status vs transfer",
                    loc="left", fontweight="bold", fontsize=9.5)

    if agree is not None and agree.shape == (2, 2):
        frac = agree.div(agree.sum(axis=1), axis=0)
        xx = np.arange(2)
        ax[1].bar(xx, frac[True].to_numpy() if True in frac else
                  frac.iloc[:, -1].to_numpy(), width=0.55,
                  color=[ORANGE, AQUA])
        ax[1].set_xticks(xx, ["STR fail", "STR pass"], fontsize=8)
        ax[1].set_ylabel("fraction that are their own best fingerprint match")
        ax[1].set_title("b  The two identity signals barely agree", loc="left",
                        fontweight="bold", fontsize=9.5)

    ax[2].bar([0, 1], S.frac_of_ceiling, width=0.55, color=[GREY, AQUA])
    for i_, r_ in enumerate(S.itertuples()):
        ax[2].text(i_, r_.frac_of_ceiling + 0.015,
                   f"{r_.frac_of_ceiling:.0%}\nn={r_.n_lines}", ha="center",
                   fontsize=7.5)
    ax[2].axhline(1.0, ls="--", color="#888", lw=1.2)
    ax[2].set_xticks([0, 1], ["all lines", "STR-passing"], fontsize=8)
    ax[2].set_ylim(0, 1.15)
    ax[2].set_ylabel("fraction of the within-lab ceiling")
    ax[2].set_title("c  Does STR filtering recover transfer?", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("STR profiling: an orthogonal, non-circular test of the "
                 "identity argument", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    dd = save_figure(fig, "str_identity", FIG,
                     source_data={"per_line": L, "selection": S},
                     script=__file__)
    print(f"figure bundle -> {dd}")


if __name__ == "__main__":
    main()
