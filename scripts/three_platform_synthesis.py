#!/usr/bin/env python3
"""Does mechanism set context-dependence? Four datasets, one ranking.

The mechanism claim has been tested piecemeal and the answers looked
inconsistent. This script puts them on one axis and asks what actually
replicates.

  Tahoe-100M    single-cell transcription, 47 lines, 379 drugs
  LINCS phase 1 bulk L1000 transcription, 70 lines, 2,834 compounds  (GSE92742)
  LINCS phase 2 bulk L1000 transcription, 30 lines, 273 compounds    (GSE70138)
  PRISM         pooled viability,        738 lines, ~1,500 compounds

The earlier Tahoe-vs-LINCS disagreement rested on phase 2 alone, where only six
mechanism classes were shared -- too few to separate "different biology" from
"too few points". Phase 1 is an order of magnitude larger and shares ten.

One confound is handled explicitly. CDI is a ratio, so a compound that does
nothing in these cells has a near-zero numerator AND denominator and lands at
0 by default, not because its response transfers. LINCS phase 1 is full of such
compounds (antihistamines, adrenergics: 35 classes sit at exactly 0.000). We
therefore report the ranking twice -- all annotated compounds, and restricted to
compounds with a detectable additive effect -- because only the second is a
statement about transferability rather than about potency.

Outputs: results/tables/three_platform_mechanism_cdi.csv
         figure bundle results/figures/13_prism/three_platform_synthesis/
"""
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
FIG = ROOT / "results" / "figures" / "13_prism"
TAB = ROOT / "results" / "tables"
MIN_N = 3
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def key(s):
    return str(s).lower().strip()


def load():
    """Per-compound CDI from each dataset, on a common column naming."""
    out = {}
    t = pd.read_csv(TAB / "drug_context_dependence.csv")
    t = t[t["moa-fine"].notna() & (t["moa-fine"] != "unclear")]
    out["Tahoe"] = pd.DataFrame({"cpd": t.drug, "moa": t["moa-fine"].map(key),
                                 "cdi": t.cdi, "active": True})
    for tag, lab in (("phase1", "LINCS-1"), ("lincs", "LINCS-2")):
        f = TAB / f"lincs_mechanism_ranking_{tag}.csv"
        if not f.exists():
            continue
        l = pd.read_csv(f)
        l = l[l.estimable & l["moa-fine"].notna() & (l["moa-fine"] != "unclear")]
        out[lab] = pd.DataFrame({"cpd": l.get("pert_iname", l.perturbation),
                                 "moa": l["moa-fine"].map(key),
                                 "cdi": l["index"],
                                 "active": l.shared > l.shared.median() * 0.25})
    f = TAB / "prism_decomposition.csv"
    if f.exists():
        p = pd.read_csv(f)
        p = p[p.moa.notna() & (p.moa.astype(str) != "nan")]
        out["PRISM"] = pd.DataFrame({"cpd": p.compound, "moa": p.moa.map(key),
                                     "cdi": p.cdi,
                                     "active": p.additive >
                                     p.additive.median() * 0.25})
    return out


def rank_table(dsets, active_only):
    cols = {}
    for name, d in dsets.items():
        dd = d[d.active] if active_only else d
        g = dd.groupby("moa").cdi.agg(["median", "size"]).query(f"size>={MIN_N}")
        cols[name] = g["median"]
    return pd.DataFrame(cols)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    dsets = load()
    print("datasets loaded: " + ", ".join(
        f"{k} ({len(v)} compounds, {v.moa.nunique()} classes)"
        for k, v in dsets.items()))

    res = {}
    for tag, act in (("all compounds", False), ("active compounds only", True)):
        R = rank_table(dsets, act)
        print(f"\n=== {tag} ===")
        names = list(R.columns)
        rows = []
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                v = R[[names[i], names[j]]].dropna()
                if len(v) < 5:
                    print(f"  {names[i]:8s} vs {names[j]:8s}  n={len(v)} "
                          f"(too few shared classes)")
                    continue
                rho = stats.spearmanr(v.iloc[:, 0], v.iloc[:, 1])
                star = "**" if rho.pvalue < 0.05 else (
                    "*" if rho.pvalue < 0.10 else "")
                print(f"  {names[i]:8s} vs {names[j]:8s}  rho={rho.statistic:+.3f}"
                      f"  p={rho.pvalue:.3f}  n={len(v)} classes {star}")
                rows.append({"a": names[i], "b": names[j], "rho": rho.statistic,
                             "p": rho.pvalue, "n": len(v), "set": tag})
        res[tag] = pd.DataFrame(rows)
        if act:
            R.to_csv(TAB / "three_platform_mechanism_cdi.csv")
            Rk = R
    A = pd.concat(res.values(), ignore_index=True)

    # is the transcriptional consensus itself consistent?
    tx = [c for c in ("Tahoe", "LINCS-1", "LINCS-2") if c in Rk.columns]
    print(f"\ntranscriptional platforms ({', '.join(tx)}) vs viability (PRISM):")
    if "PRISM" in Rk.columns and len(tx) >= 2:
        z = Rk[tx].apply(lambda c: c.rank(pct=True))
        cons = z.mean(axis=1, skipna=True)
        v = pd.DataFrame({"transcription": cons,
                          "viability": Rk.PRISM}).dropna()
        if len(v) >= 5:
            rho = stats.spearmanr(v.transcription, v.viability)
            print(f"  consensus transcription vs PRISM viability: "
                  f"rho={rho.statistic:+.3f}, p={rho.pvalue:.3f}, n={len(v)}")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8), constrained_layout=True)

    lab = [f"{r.a}\nvs {r.b}" for r in A[A.set == "active compounds only"]
           .itertuples()]
    val = A[A.set == "active compounds only"].rho.to_numpy()
    nn = A[A.set == "active compounds only"].n.to_numpy()
    col = [AQUA if v > 0 else ORANGE for v in val]
    ax[0].bar(np.arange(len(val)), val, color=col, width=0.6)
    ax[0].axhline(0, color="#444444", lw=0.9)
    ax[0].set_xticks(np.arange(len(val)), lab, fontsize=6.5, rotation=45,
                     ha="right")
    for i, (v, k) in enumerate(zip(val, nn)):
        ax[0].text(i, v + (0.03 if v > 0 else -0.06), f"n={k}", ha="center",
                   fontsize=6)
    ax[0].set_ylabel("Spearman rho of mechanism CDI")
    ax[0].set_title("A  Cross-dataset agreement (active compounds)", loc="left",
                    fontweight="bold", fontsize=10)

    if "LINCS-1" in Rk.columns:
        v = Rk[["Tahoe", "LINCS-1"]].dropna()
        rho = stats.spearmanr(v.Tahoe, v["LINCS-1"])
        ax[1].scatter(v.Tahoe, v["LINCS-1"], s=40, color=BLUE, edgecolors="none")
        for m, r in v.iterrows():
            ax[1].annotate(m[:22], (r.Tahoe, r["LINCS-1"]), fontsize=6,
                           xytext=(3, 2), textcoords="offset points")
        ax[1].set_xlabel("mechanism CDI — Tahoe-100M (single-cell)")
        ax[1].set_ylabel("mechanism CDI — LINCS phase 1 (bulk L1000)")
        ax[1].set_title(f"B  Transcription replicates\nrho={rho.statistic:+.2f}, "
                        f"n={len(v)}", loc="left", fontweight="bold", fontsize=10)

    if "PRISM" in Rk.columns:
        v = Rk[["Tahoe", "PRISM"]].dropna()
        rho = stats.spearmanr(v.Tahoe, v.PRISM)
        ax[2].scatter(v.Tahoe, v.PRISM, s=40, color=VIOLET, edgecolors="none")
        for m, r in v.iterrows():
            ax[2].annotate(m[:22], (r.Tahoe, r.PRISM), fontsize=6,
                           xytext=(3, 2), textcoords="offset points")
        ax[2].set_xlabel("mechanism CDI — Tahoe-100M (transcription)")
        ax[2].set_ylabel("mechanism CDI — PRISM (viability)")
        ax[2].set_title(f"C  Does it cross the readout gap?\n"
                        f"rho={rho.statistic:+.2f}, n={len(v)}", loc="left",
                        fontweight="bold", fontsize=10)
    fig.suptitle("Does drug mechanism set context-dependence — in which "
                 "datasets, and across which readouts?", fontsize=11, x=0.01,
                 ha="left")
    d = save_figure(fig, "three_platform_synthesis", FIG,
                    source_data={"mechanism_cdi": Rk.reset_index(),
                                 "agreements": A}, script=__file__)
    print(f"\nfigure bundle -> {d}")


if __name__ == "__main__":
    main()
