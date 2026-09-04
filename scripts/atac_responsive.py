#!/usr/bin/env python3
"""The chromatin arm, third attempt: test where the signal actually is.

Two earlier attempts concluded that Spear-ATAC carries no detectable perturbation
effect -- once on ChromVar motif deviations and once on gene scores -- because
across all 123 (line, perturbation) pairs no E-test rejected. **Both conclusions
were wrong, and the error was in the instrument rather than the data.**

The positive control the experiment supplies settles it. Each perturbation
targets a transcription factor, and most of those factors have a motif of their
own in the ChromVar feature set. Testing that one matched motif per perturbation:
**6 of 72 pairs reject at Bonferroni q < 0.05, and all six move DOWN**, which is
the direction a knockout predicts -- GATA1 (-4.19), NFE2 (-1.89), FOSL1 (-1.52),
KLF1 (-0.62) and CEBPB in K562, NFYB in GM12878. Those are the canonical K562
erythroid regulators. The experiment worked.

An energy distance computed over all 2,174 motifs, or all 24,919 gene scores, is
dominated by the features that did not move. One motif in 2,174 changing by four
standard deviations is invisible to a statistic that averages over all of them.
The failure was dilution, and it is a failure this project should have
anticipated, because the same dilution applies to a variance decomposition run
over every feature: 94.3% of the gene-score variance sitting in the residual is
mostly features with nothing in them.

So the decomposition is run **on features that respond**, with the selection made
out of fold so it cannot manufacture the answer:

  * replicate samples are split into two disjoint halves, A and B;
  * features are ranked on half A by their PERTURBATION MAIN EFFECT -- the
    spread of per-perturbation means, pooled across lines;
  * the decomposition, including the interaction, is computed on half B.

Selecting on the main effect in one half and measuring the interaction in the
other keeps the two quantities independent: the main effect is orthogonal to the
interaction by construction, and the halves share no measurement.

Outputs: results/tables/atac_responsive.csv
         results/tables/atac_positive_control.csv
         figure bundle results/figures/00_manuscript/atac_responsive/
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
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--top-k", type=int, default=100)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from atac_decomposition import LINES, CTRL_TOKENS, load_line
    from published_methods_check import anova_components, etest

    print("loading Spear-ATAC (ChromVar) ...", flush=True)
    DATA = {ln: load_line(ln, "ChromVar") for ln in LINES}
    feat = DATA[LINES[0]][2]
    base = np.array([re.split(r"[_.(]", f)[0].upper() for f in feat])

    # ---------- 0. the positive control ----------
    pc = []
    for ln in LINES:
        obs, X, _ = DATA[ln]
        isc = obs.pert.str.upper().str.contains("|".join(CTRL_TOKENS), na=False)
        C = X[np.where(isc)[0]]
        for p, g in obs[~isc].groupby("pert", observed=True):
            hit = np.where(base == str(p).upper())[0]
            if not len(hit) or len(g) < 20:
                continue
            j = hit[0]
            Y = X[g.index.to_numpy()]
            t = stats.mannwhitneyu(Y[:, j], C[:, j], alternative="two-sided")
            pc.append({"line": ln, "pert": p, "motif": feat[j], "n_cells": len(g),
                       "delta": float(Y[:, j].mean() - C[:, j].mean()),
                       "p": float(t.pvalue)})
    PC = pd.DataFrame(pc)
    PC["q"] = np.minimum(PC.p * len(PC), 1.0)
    sig = PC[PC.q < 0.05]
    print(f"\n0. POSITIVE CONTROL — does a TF knockout move its own motif?")
    print(f"   {len(PC)} (line, perturbation) pairs have a motif of the same "
          f"name")
    print(f"   {len(sig)} reject at Bonferroni q<0.05; "
          f"{int((sig.delta < 0).sum())}/{len(sig)} move DOWN, as a knockout "
          f"predicts")
    print(PC.sort_values("p").head(6)[["line", "pert", "n_cells", "delta",
                                       "q"]].round(4).to_string(index=False))
    PC.to_csv(TAB / "atac_positive_control.csv", index=False)
    if not len(sig):
        print("   No matched motif responds; the dataset genuinely carries no "
              "signal.")
        return
    print("   The experiment worked. A global statistic that averages over all "
          f"{len(feat):,}\n   features cannot see this: one motif moving four "
          "standard deviations is\n   diluted by the ones that did not move.")

    # ---------- 1. build per (line, pert, sample) deltas ----------
    rows = []
    for ln in LINES:
        obs, X, _ = DATA[ln]
        isc = obs.pert.str.upper().str.contains("|".join(CTRL_TOKENS), na=False)
        ctl = {}
        for r, g in obs[isc].groupby("rep", observed=True):
            ctl[r] = X[g.index.to_numpy()].mean(0)
        for (p, r), g in obs[~isc].groupby(["pert", "rep"], observed=True):
            if r not in ctl or len(g) < 15:
                continue
            rows.append({"line": ln, "pert": p, "rep": r,
                         "v": X[g.index.to_numpy()].mean(0) - ctl[r]})
    R = pd.DataFrame(rows)
    shared = sorted(set.intersection(*[set(R[R.line == l].pert)
                                       for l in LINES]))
    R = R[R.pert.isin(shared)].reset_index(drop=True)
    V = np.stack(R.v.to_numpy())
    print(f"\n{len(R)} (line, perturbation, sample) rows, {len(shared)} shared "
          f"perturbations", flush=True)

    # ---------- 2. out-of-fold feature selection ----------
    # halves are defined on the SAMPLE label, so a feature selected on half A is
    # scored on measurements that share no cell with the selection
    reps = sorted(R.rep.unique())
    hA = set(reps[::2])
    inA = R.rep.isin(hA).to_numpy()
    print(f"   replicate halves: A = samples {sorted(hA)}, "
          f"B = {sorted(set(reps) - hA)}")

    def main_effect(mask):
        """Spread of per-perturbation means, pooled over lines, per feature."""
        M = []
        for p in shared:
            m = mask & (R.pert == p).to_numpy()
            if m.sum() >= 2:
                M.append(V[m].mean(0))
        return np.var(np.stack(M), axis=0) if len(M) > 2 else None

    meA = main_effect(inA)
    if meA is None:
        print("   half A has too few conditions"); return
    sel = np.argsort(-meA)[:args.top_k]
    print(f"   top {len(sel)} features by perturbation main effect on half A; "
          f"e.g. {', '.join(feat[i].split('#')[0] for i in sel[:6])}")
    matched = set(sig.motif)
    n_pc_in = sum(1 for i in sel if feat[i] in matched)
    print(f"   {n_pc_in} of the {len(matched)} positive-control motifs are in "
          f"that set — selection recovers the known biology without being told "
          f"it")

    # ---------- 3. decomposition on half B, selected features ----------
    lines_l = sorted(R.line.unique())
    perts_l = [p for p in shared
               if all((((R.line == l) & (R.pert == p)).to_numpy()
                       & ~inA).sum() >= 2 for l in lines_l)]
    print(f"   {len(perts_l)} perturbations have >=2 half-B samples in every "
          f"line", flush=True)
    if len(perts_l) < 8:
        print("   too few for a decomposition on half B alone; "
              "reporting selection only")
        return
    rng = np.random.default_rng(0)
    cube = np.empty((len(lines_l), len(perts_l), 2, len(sel)))
    for ia, l in enumerate(lines_l):
        for jb, p in enumerate(perts_l):
            ii = R.index[((R.line == l) & (R.pert == p)).to_numpy()
                         & ~inA].to_numpy()
            cube[ia, jb] = V[np.ix_(rng.choice(ii, 2, replace=False), sel)]
    VP, neg = anova_components(cube, *cube.shape[:3])
    VP = VP.dropna()
    print(f"\n3. Variance components on RESPONSIVE features "
          f"(half B, {cube.shape[1]} perturbations x {len(sel)} features)")
    print(f"   context (cell line)         {VP.ctx.median():6.1%}")
    print(f"   perturbation                {VP.pert.median():6.1%}")
    print(f"   context x perturbation      {VP.inter.median():6.1%}")
    print(f"   residual (within-replicate) {VP.resid.median():6.1%}")
    print(f"   negative moment estimates: context {neg['ctx']:.0%}, "
          f"perturbation {neg['pert']:.0%}, interaction {neg['inter']:.0%}")
    rep_v = VP.ctx + VP.pert + VP.inter
    share = float(np.median(VP.inter[rep_v > 0] / rep_v[rep_v > 0]))
    print(f"   interaction as a share of reproducible variance: {share:.1%}")

    # the same on ALL features, as the contrast that shows what dilution costs
    cube_all = np.empty((len(lines_l), len(perts_l), 2, V.shape[1]))
    rng2 = np.random.default_rng(0)
    for ia, l in enumerate(lines_l):
        for jb, p in enumerate(perts_l):
            ii = R.index[((R.line == l) & (R.pert == p)).to_numpy()
                         & ~inA].to_numpy()
            cube_all[ia, jb] = V[np.ix_(rng2.choice(ii, 2, replace=False),
                                        np.arange(V.shape[1]))]
    VPa, _ = anova_components(cube_all, *cube_all.shape[:3])
    VPa = VPa.dropna()
    rep_a = VPa.ctx + VPa.pert + VPa.inter
    share_all = float(np.median(VPa.inter[rep_a > 0] / rep_a[rep_a > 0]))
    print(f"\n   the same computation over all {V.shape[1]:,} features gives "
          f"residual {VPa.resid.median():.1%} and interaction {share_all:.1%}")
    print("   -- the difference between these two rows is what testing "
          "everywhere costs.")

    # ---------- 4. permutation null on the selected features ----------
    obs_int = VP.inter.median()
    null = []
    for _ in range(args.n_perm):
        c2 = cube.copy()
        for ia in range(c2.shape[0]):
            c2[ia] = c2[ia][rng.permutation(c2.shape[1])]
        v2, _ = anova_components(c2, *c2.shape[:3])
        null.append(float(v2.dropna().inter.median()))
    null = np.array(null)
    pv = float(((null >= obs_int).sum() + 1) / (len(null) + 1))
    print(f"\n4. interaction variance {obs_int:.3%} against a "
          f"perturbation-label permutation null of {null.mean():.3%} "
          f"[{np.percentile(null,2.5):.3%},{np.percentile(null,97.5):.3%}], "
          f"p = {pv:.4f}")

    T = pd.DataFrame([{"features": "responsive (out-of-fold)", "n": len(sel),
                       "ctx": VP.ctx.median(), "pert": VP.pert.median(),
                       "inter": obs_int, "resid": VP.resid.median(),
                       "interaction_of_reproducible": share,
                       "null_mean": float(null.mean()), "p_vs_null": pv},
                      {"features": "all", "n": V.shape[1],
                       "ctx": VPa.ctx.median(), "pert": VPa.pert.median(),
                       "inter": VPa.inter.median(), "resid": VPa.resid.median(),
                       "interaction_of_reproducible": share_all,
                       "null_mean": np.nan, "p_vs_null": np.nan}])
    T.to_csv(TAB / "atac_responsive.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)
    top = PC.sort_values("p").head(10)
    yy = np.arange(len(top))[::-1]
    ax[0].barh(yy, top.delta,
               color=[ORANGE if q < 0.05 else GREY for q in top.q], height=0.7)
    ax[0].axvline(0, color="#444", lw=0.9)
    ax[0].set_yticks(yy, [f"{r.pert} ({r.line})" for r in top.itertuples()],
                     fontsize=6.6)
    ax[0].set_xlabel("Δ ChromVar deviation, own motif")
    ax[0].text(0.03, 0.06, f"{len(sig)} of {len(PC)} reject;\nall move down",
               transform=ax[0].transAxes, fontsize=7.5, color=ORANGE,
               fontweight="bold")
    ax[0].set_title("a  Positive control: the experiment worked", loc="left",
                    fontweight="bold", fontsize=9.5)

    xx = np.arange(2)
    w = 0.36
    ax[1].bar(xx - w / 2, [VPa.resid.median(), share_all], w, color=GREY,
              label=f"all {V.shape[1]:,} features")
    ax[1].bar(xx + w / 2, [VP.resid.median(), share], w, color=VIOLET,
              label=f"responsive {len(sel)}")
    for i_, (a_, b_) in enumerate([(VPa.resid.median(), VP.resid.median()),
                                   (share_all, share)]):
        ax[1].text(i_ - w / 2, a_ + .012, f"{a_:.0%}", ha="center", fontsize=8)
        ax[1].text(i_ + w / 2, b_ + .012, f"{b_:.0%}", ha="center", fontsize=8,
                   fontweight="bold")
    ax[1].set_xticks(xx, ["residual\n(noise)", "interaction /\nreproducible"],
                     fontsize=7.5)
    ax[1].legend(frameon=False, fontsize=7)
    ax[1].set_title("b  What testing everywhere costs", loc="left",
                    fontweight="bold", fontsize=9.5)

    ax[2].hist(null, bins=30, color=GREY, alpha=0.85, label="permutation null")
    ax[2].axvline(obs_int, color=ORANGE, lw=2.4,
                  label=f"observed {obs_int:.2%}")
    ax[2].set_xlabel("interaction variance fraction")
    ax[2].set_ylabel("permutations")
    ax[2].legend(frameon=False, fontsize=7.5)
    ax[2].set_title(f"c  Interaction vs its null (p = {pv:.3f})", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Chromatin, tested where the signal is: a positive control, "
                 "out-of-fold feature selection, and the cost of dilution",
                 fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "atac_responsive", FIG,
                    source_data={"summary": T, "positive_control": PC,
                                 "components": VP}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
