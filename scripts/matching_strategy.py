#!/usr/bin/env python3
"""What does the cell-line matching strategy cost you?

Every cross-dataset comparison in this field begins by deciding that a line in
one atlas is "the same" as a line in another. That decision is almost always made
from an identifier -- a name, a DepMap ID, a COSMIC ID -- and then never checked.
Ben-David et al. (2018) showed the assumption can fail: laboratories culture
divergent strains, and a median 19% of non-silent mutations differ between CCLE
and GDSC for nominally identical lines.

This script measures what the choice is worth, by running the same cross-lab
comparison (PRISM at the Broad vs GDSC at the Sanger) under a ladder of matching
strategies of increasing stringency:

  random          lines paired at random — the null, and the floor
  tissue          paired at random WITHIN the same primary tissue — tests how
                  much apparent agreement is lineage rather than identity
  identifier      the curated COSMIC->DepMap match, i.e. standard practice
  fingerprint     the identifier match, kept only where the line is its own
                  reciprocal best hit by drug-response fingerprint
  best-hit        each line paired with its most similar partner regardless of
                  identifier — an upper bound, and deliberately circular, shown
                  only to bound how much room the data allows

The tissue rung matters most for interpretation. If identifier-matched pairs
barely beat same-tissue random pairs, then what looks like "the same cell line
behaves the same way" is really "cells from the same lineage behave the same
way", and cross-atlas identity is carrying almost no information.

A baseline-expression rung is deliberately absent, and the reason is worth
stating: LINCS Level 4 is z-scored within plate, which removes the cell-line
baseline that expression matching would key on, and neither PRISM nor GDSC ships
expression at all. The response fingerprint is what remains measurable in every
dataset here, and it is arguably the more relevant quantity anyway -- it asks
whether two cultures behave alike pharmacologically, which is the property a
cross-atlas analysis actually relies on.

Outputs: results/tables/matching_strategy.csv
         figure bundle results/figures/16_crosslab/matching_strategy/
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
GD = ROOT / "data" / "external" / "gdsc"
FIG = ROOT / "results" / "figures" / "16_crosslab"
TAB = ROOT / "results" / "tables"
MIN_SHARED = 20
MIN_CPD_FP = 8
N_RAND = 20
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def build():
    from importlib import import_module
    m = import_module("cross_lab_reproducibility")
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-rand", type=int, default=N_RAND)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    m = build()

    print("building residuals ...", flush=True)
    p_full, p_h1, p_h2 = m.prism_residuals()
    G = m.gdsc_residuals()
    smp = pd.read_csv(PR / "sample_info.tsv", sep="\t", low_memory=False)
    dep = smp.cross_references.astype(str).str.extract(r"DepMap;\s*(ACH-\d+)")[0]
    cos = {}
    for d, xr in zip(dep, smp.cross_references.astype(str)):
        if not isinstance(d, str):
            continue
        for cid in re.findall(r"Cosmic(?:-CLP)?;\s*(\d+)", xr):
            cos.setdefault(int(cid), d)
    gb = {}
    for tag in ("GDSC1", "GDSC2"):
        g = G[tag].copy()
        g["dep"] = g.COSMIC_ID.map(cos)
        g = g[g.dep.notna()]
        for k, gg in g.groupby("k", observed=True):
            v = gg.groupby("dep").Z_SCORE.mean()
            gb.setdefault(k, []).append(v)
    gboth = {k: pd.concat(v).groupby(level=0).mean() for k, v in gb.items()}
    P = {norm(k): v for k, v in p_full.items()}
    shared_cpd = sorted(set(P) & set(gboth))
    print(f"  {len(shared_cpd)} shared compounds", flush=True)

    # per-line fingerprints over the shared compounds
    def fingerprints(store):
        F = {}
        for c in shared_cpd:
            s = store.get(c)
            if s is None:
                continue
            for ln, val in s.items():
                F.setdefault(ln, {})[c] = float(val)
        return {k: v for k, v in F.items() if len(v) >= MIN_CPD_FP}

    FA, FB = fingerprints(P), fingerprints(gboth)
    both = sorted(set(FA) & set(FB))
    print(f"  {len(FA)} PRISM lines, {len(FB)} GDSC lines, "
          f"{len(both)} identifier-matched", flush=True)

    def fp_corr(x, y):
        common = sorted(set(FA[x]) & set(FB[y]))
        if len(common) < MIN_CPD_FP:
            return np.nan
        u = np.array([FA[x][c] for c in common])
        v = np.array([FB[y][c] for c in common])
        if np.std(u) == 0 or np.std(v) == 0:
            return np.nan
        return float(stats.spearmanr(u, v).statistic)

    print("computing the similarity matrix ...", flush=True)
    lb = sorted(FB)
    Mx = pd.DataFrame(np.nan, index=both, columns=lb)
    for x in both:
        for y in lb:
            Mx.loc[x, y] = fp_corr(x, y)

    ci = pd.read_csv(PR / "secondary-screen-cell-line-info.csv")
    tis = dict(zip(ci.depmap_id, ci.primary_tissue.astype(str)))

    rows = []
    # identifier
    for x in both:
        r = Mx.loc[x, x]
        if np.isfinite(r):
            rows.append({"strategy": "identifier\n(standard practice)",
                         "line": x, "r": r})
    # fingerprint-validated subset of the identifier match
    validated = []
    for x in both:
        s = Mx.loc[x].dropna()
        if len(s) and s.idxmax() == x:
            validated.append(x)
            rows.append({"strategy": "identifier +\nfingerprint-validated",
                         "line": x, "r": float(s.max())})
    # best-hit (upper bound; circular by construction)
    for x in both:
        s = Mx.loc[x].dropna()
        if len(s):
            rows.append({"strategy": "best hit\n(upper bound)", "line": x,
                         "r": float(s.max())})
    # random and tissue-matched nulls
    for _ in range(args.n_rand):
        perm = rng.permutation(len(both))
        for i, x in enumerate(both):
            y = both[perm[i]]
            if y == x:
                continue
            r = Mx.loc[x, y]
            if np.isfinite(r):
                rows.append({"strategy": "random pairing", "line": x, "r": r})
        for x in both:
            same = [y for y in both if y != x
                    and tis.get(y, "?") == tis.get(x, "?") != "?"]
            if not same:
                continue
            y = same[rng.integers(0, len(same))]
            r = Mx.loc[x, y]
            if np.isfinite(r):
                rows.append({"strategy": "same tissue,\nrandom line", "line": x,
                             "r": r})
    R = pd.DataFrame(rows)
    R.to_csv(TAB / "matching_strategy.csv", index=False)

    order = ["random pairing", "same tissue,\nrandom line",
             "identifier\n(standard practice)",
             "identifier +\nfingerprint-validated", "best hit\n(upper bound)"]
    order = [o for o in order if (R.strategy == o).any()]
    S = R.groupby("strategy").r.agg(["median", "mean", "size"]).reindex(order)
    print("\n=== fingerprint similarity by matching strategy ===")
    print(S.round(3).to_string())
    idm = R[R.strategy == order[2]].r
    tism = R[R.strategy == order[1]].r
    if len(idm) and len(tism):
        u = stats.mannwhitneyu(idm, tism)
        print(f"\nidentifier vs same-tissue-random: {idm.median():.3f} vs "
              f"{tism.median():.3f}, p={u.pvalue:.2e}")
        print("  The gap between these two is what cross-atlas cell-line "
              "IDENTITY buys\n  over merely matching lineage. If it is small, "
              "identity is carrying little.")
    print(f"\n{len(validated)}/{len(both)} identifier matches survive "
          f"fingerprint validation ({len(validated)/max(len(both),1):.0%})")

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8), constrained_layout=True)
    data = [R[R.strategy == o].r.dropna() for o in order]
    cols = ["#bdbdbd", AQUA, ORANGE, VIOLET, BLUE][:len(order)]
    bp = ax[0].boxplot(data, showfliers=False, patch_artist=True,
                       medianprops=dict(color="black", lw=1.6))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.8)
    ax[0].set_xticks(range(1, len(order) + 1), order, fontsize=6.5)
    ax[0].axhline(0, color="#444444", lw=0.9)
    for i, d_ in enumerate(data):
        if len(d_):
            ax[0].text(i + 1, np.median(d_) + 0.02, f"{np.median(d_):.2f}",
                       ha="center", fontsize=8.5, fontweight="bold")
    ax[0].set_ylabel("response-fingerprint similarity (Spearman r)")
    ax[0].set_title("A  What each matching strategy buys", loc="left",
                    fontweight="bold", fontsize=10)

    for o, c in zip(order, cols):
        d_ = R[R.strategy == o].r.dropna()
        if len(d_) > 5:
            ax[1].hist(d_, bins=40, histtype="step", lw=2, color=c,
                       density=True, label=o.replace("\n", " "))
    ax[1].axvline(0, color="#444444", lw=0.9)
    ax[1].set_xlabel("fingerprint similarity"); ax[1].set_ylabel("density")
    ax[1].legend(frameon=False, fontsize=6.5)
    ax[1].set_title("B  Distributions", loc="left", fontweight="bold",
                    fontsize=10)

    idr = R[R.strategy == order[2]].set_index("line").r
    keep = idr.index.isin(validated)
    ax[2].scatter(np.where(keep, 1, 0) + rng.normal(0, 0.06, len(idr)),
                  idr.to_numpy(), s=14, alpha=0.45,
                  color=[VIOLET if k else ORANGE for k in keep],
                  edgecolors="none")
    ax[2].set_xticks([0, 1], [f"failed validation\n(n={int((~keep).sum())})",
                              f"validated\n(n={int(keep.sum())})"], fontsize=8)
    ax[2].axhline(0, color="#444444", lw=0.9)
    ax[2].set_ylabel("fingerprint similarity of the identifier match")
    ax[2].set_title(f"C  {len(validated)}/{len(both)} identifier matches "
                    f"survive\n({len(validated)/max(len(both),1):.0%})",
                    loc="left", fontweight="bold", fontsize=10)
    fig.suptitle("Matching cell lines across atlases: what identifiers assume, "
                 "and what the data supports", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "matching_strategy", FIG,
                    source_data={"per_pair": R, "summary": S.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
