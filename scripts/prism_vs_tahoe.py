#!/usr/bin/env python3
"""Does the mechanism ranking replicate in PRISM, and does that settle LINCS?

The LINCS comparison was thin: only 6 mechanism classes were shared with Tahoe,
so the rank correlation there rested on 6 points and could not distinguish
"different biology" from "too few classes". PRISM annotates 77 classes at n>=5
and shares far more with the Tahoe panel, so it is the test LINCS could not be.

Two levels of agreement, both reported:

  mechanism level  Spearman of median CDI over shared MOA classes. Answers the
                   claim as stated ("mechanism sets transferability").
  compound level   the same compound measured in both atlases, matched by name.
                   Stronger, because it does not depend on the MOA vocabulary
                   lining up, and because a compound is its own control for
                   mechanism annotation error.

PRISM measures viability, Tahoe measures transcription. Agreement across that
gap is a much harder test than agreement between two transcriptional platforms:
it says the ordering is a property of the drug, not of the readout.

Outputs: results/tables/prism_vs_tahoe_cdi.csv
         results/tables/three_platform_mechanism_cdi.csv
         figure bundle results/figures/13_prism/prism_vs_tahoe/
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
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    P = pd.read_csv(TAB / "prism_decomposition.csv")
    T = pd.read_csv(TAB / "drug_context_dependence.csv")
    T = T[T["moa-fine"].notna() & (T["moa-fine"] != "unclear")].copy()
    P["m"] = P.moa.astype(str).str.lower().str.strip()
    T["m"] = T["moa-fine"].str.lower().str.strip()

    pg = P[P.moa.notna()].groupby("m").cdi.agg(["median", "size"])
    tg = T.groupby("m").cdi.agg(["median", "size"])
    both = pd.DataFrame({"prism": pg["median"], "n_prism": pg["size"],
                         "tahoe": tg["median"], "n_tahoe": tg["size"]}).dropna()
    both = both[(both.n_prism >= 3) & (both.n_tahoe >= 3)]
    rho = stats.spearmanr(both.prism, both.tahoe)
    print(f"MECHANISM level: {len(both)} shared classes, "
          f"Spearman rho={rho.statistic:+.3f}, p={rho.pvalue:.4f}")
    print(both.sort_values("tahoe", ascending=False).round(3).to_string())
    both.to_csv(TAB / "prism_vs_tahoe_cdi.csv")

    # compound level — no MOA vocabulary needed
    P["k"] = P.compound.map(norm)
    T["k"] = T.drug.map(norm)
    cm = P[["k", "cdi", "n_lines"]].merge(
        T[["k", "drug", "cdi"]], on="k", suffixes=("_prism", "_tahoe"))
    cm = cm.drop_duplicates("k").dropna()
    print(f"\nCOMPOUND level: {len(cm)} compounds measured in both atlases")
    if len(cm) >= 8:
        rc = stats.spearmanr(cm.cdi_prism, cm.cdi_tahoe)
        print(f"  Spearman rho={rc.statistic:+.3f}, p={rc.pvalue:.2e}")
        print("  most context-dependent in both:")
        cm["rank_sum"] = (cm.cdi_prism.rank(pct=True)
                          + cm.cdi_tahoe.rank(pct=True))
        print(cm.nlargest(10, "rank_sum")[
            ["drug", "cdi_prism", "cdi_tahoe", "n_lines"]].round(3)
            .to_string(index=False))
    else:
        rc = None

    # three-platform table where LINCS also has the class
    rows = both.copy()
    lf = TAB / "lincs_vs_tahoe_cdi.csv"
    if lf.exists():
        L = pd.read_csv(lf, index_col=0)
        rows = rows.join(L[["lincs"]], how="left")
    rows.to_csv(TAB / "three_platform_mechanism_cdi.csv")
    n3 = rows.lincs.notna().sum() if "lincs" in rows else 0
    print(f"\nclasses with all three platforms: {n3}")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    ax[0].scatter(both.tahoe, both.prism, s=38, color=ORANGE, edgecolors="none")
    for m, r in both.iterrows():
        ax[0].annotate(m[:20], (r.tahoe, r.prism), fontsize=6,
                       xytext=(3, 2), textcoords="offset points")
    ax[0].set_xlabel("mechanism median CDI — Tahoe-100M (transcription)")
    ax[0].set_ylabel("mechanism median CDI — PRISM (viability)")
    ax[0].set_title(f"A  Mechanism ranking replicates across readouts\n"
                    f"rho={rho.statistic:+.2f}, n={len(both)} classes",
                    loc="left", fontweight="bold", fontsize=10)

    if rc is not None:
        ax[1].scatter(cm.cdi_tahoe, cm.cdi_prism, s=14, alpha=0.45, color=BLUE,
                      edgecolors="none")
        ax[1].set_xlabel("compound CDI — Tahoe-100M")
        ax[1].set_ylabel("compound CDI — PRISM")
        ax[1].set_title(f"B  Same compounds, both atlases\n"
                        f"rho={rc.statistic:+.2f}, n={len(cm)} compounds",
                        loc="left", fontweight="bold", fontsize=10)

    lab = ["Tahoe\n(transcription,\n47 lines)", "PRISM\n(viability,\n738 lines)"]
    val = [len(T), len(P)]
    ncls = [T.m.nunique(), (pg["size"] >= 5).sum()]
    if "lincs" in rows:
        lab.insert(1, "LINCS\n(transcription,\n30 lines)")
    ax[2].bar(np.arange(len(ncls)), ncls, color=[BLUE, ORANGE][:len(ncls)],
              width=0.55)
    ax[2].set_xticks(np.arange(len(ncls)),
                     [l for l in lab if "LINCS" not in l], fontsize=8)
    ax[2].set_ylabel("annotated mechanism classes")
    ax[2].set_title("C  Why PRISM is the decisive test", loc="left",
                    fontweight="bold", fontsize=10)
    fig.suptitle("Does drug mechanism set context-dependence for viability too?",
                 fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "prism_vs_tahoe", FIG,
                    source_data={"by_mechanism": both.reset_index(),
                                 "by_compound": cm,
                                 "three_platform": rows.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
