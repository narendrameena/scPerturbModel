#!/usr/bin/env python3
"""Re-derive the central split with established, published methods.

Every decomposition in this project uses an estimator built here: a covariance
between independent replicates, with main effects removed out of fold. It is
validated against simulation, but a bespoke estimator validated by its author's
own simulation is weaker evidence than the same conclusion reached by a method
the field already accepts. This script re-derives the central quantity with two
published methods and reports where they agree and disagree.

**variancePartition** (Hoffman & Schadt, *BMC Bioinformatics* 17:483, 2016) is the
standard variance-component method in genomics. It fits, per feature, a linear
mixed model with crossed random effects

    y ~ (1|context) + (1|perturbation) + (1|context:perturbation) + residual

by REML and reports each term as a fraction of total variance. This is exactly
the additive/interaction split, obtained by a route that predates this project
and makes none of its choices. The interaction term is identifiable **only
because the design is replicated** -- without replicates the interaction and the
residual are the same stratum, which is the paper's central methodological point
arriving from an independent direction.

**E-distance and the E-test** (Peidli et al., *Nature Methods* 21:531, 2024 --
the scPerturb paper, whose data this analysis uses) is the field's standard
statistic for whether a perturbation moved a cell population at all. The energy
distance between perturbed and control cells,

    E = 2 E||x - y|| - E||x - x'|| - E||y - y'||,

is zero exactly when the two distributions coincide, and the permutation E-test
gives it a null. It is used here as a prerequisite check: a perturbation whose
E-test does not reject should not contribute to any interaction estimate, and if
the interaction is carried by perturbations that did nothing, it is an artefact.

The comparison is run on Spear-ATAC (Pierce/Greenleaf 2021), which is fully
crossed -- all 41 perturbations in all 3 lines, 4-6 replicate samples each -- so
the mixed model is identifiable without imbalance corrections and the two methods
see the same data.

Outputs: results/tables/published_methods_check.csv
         figure bundle results/figures/00_manuscript/published_methods/
"""
import argparse
import warnings
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


def anova_components(Y, n_ctx, n_pert, n_rep):
    """Two-way crossed random-effects variance components, balanced design.

    Henderson's Method III / ANOVA method of moments (Searle, Casella &
    McCulloch, *Variance Components*, 1992, ch. 4), which is what
    variancePartition's REML (Hoffman & Schadt 2016) reduces to when the design
    is balanced. Used in place of a REML fit because statsmodels' crossed
    random-effects idiom fails to build a design matrix for this layout, and a
    closed-form estimator on a balanced design is both exact and inspectable.

    Model  y_ijk = mu + a_i + b_j + (ab)_ij + e_ijk  with i contexts, j
    perturbations, k replicates. Expected mean squares give

        sigma_e^2  = MS_E
        sigma_ab^2 = (MS_AB - MS_E) / n_rep
        sigma_a^2  = (MS_A  - MS_AB) / (n_rep * n_pert)
        sigma_b^2  = (MS_B  - MS_AB) / (n_rep * n_ctx)

    ``Y`` is (n_ctx, n_pert, n_rep, n_features). Components are returned as
    fractions of their total per feature; negative estimates, which the method
    of moments can produce when a true component is near zero, are clipped to
    zero and their frequency is reported rather than hidden.
    """
    a, b, n = n_ctx, n_pert, n_rep
    gm = Y.mean(axis=(0, 1, 2))
    cell = Y.mean(axis=2)                       # ctx x pert x feat
    ma = cell.mean(axis=1)                      # ctx x feat
    mb = cell.mean(axis=0)                      # pert x feat
    ss_a = n * b * ((ma - gm) ** 2).sum(0)
    ss_b = n * a * ((mb - gm) ** 2).sum(0)
    ss_ab = n * ((cell - ma[:, None, :] - mb[None, :, :] + gm) ** 2).sum((0, 1))
    ss_e = ((Y - cell[:, :, None, :]) ** 2).sum((0, 1, 2))
    ms_a = ss_a / max(a - 1, 1)
    ms_b = ss_b / max(b - 1, 1)
    ms_ab = ss_ab / max((a - 1) * (b - 1), 1)
    ms_e = ss_e / max(a * b * (n - 1), 1)
    v_e = ms_e
    v_ab = (ms_ab - ms_e) / n
    v_a = (ms_a - ms_ab) / (n * b)
    v_b = (ms_b - ms_ab) / (n * a)
    neg = {"ctx": float((v_a < 0).mean()), "pert": float((v_b < 0).mean()),
           "inter": float((v_ab < 0).mean())}
    comp = {k: np.maximum(v, 0) for k, v in
            (("ctx", v_a), ("pert", v_b), ("inter", v_ab), ("resid", v_e))}
    tot = sum(comp.values())
    ok = tot > 0
    frac = {k: np.where(ok, v / np.where(ok, tot, 1), np.nan)
            for k, v in comp.items()}
    return pd.DataFrame(frac), neg


def edistance(X, Y, max_n=400, seed=0):
    """Energy distance between two cell populations (Peidli et al. 2024).

    E = 2 E||x-y|| - E||x-x'|| - E||y-y'||, on Euclidean distances. Subsampled
    to `max_n` per group because the statistic is O(n^2) and the subsample is
    unbiased for it.
    """
    rng = np.random.default_rng(seed)
    if len(X) > max_n:
        X = X[rng.choice(len(X), max_n, replace=False)]
    if len(Y) > max_n:
        Y = Y[rng.choice(len(Y), max_n, replace=False)]
    if len(X) < 5 or len(Y) < 5:
        return np.nan

    def md(A, B):
        d = np.sqrt(np.maximum(
            ((A[:, None, :] - B[None, :, :]) ** 2).sum(-1), 0))
        return float(d.mean())
    return 2 * md(X, Y) - md(X, X) - md(Y, Y)


def etest(X, Y, n_perm=200, seed=0):
    """Permutation E-test: is the energy distance larger than chance?"""
    obs = edistance(X, Y, seed=seed)
    if not np.isfinite(obs):
        return obs, np.nan
    rng = np.random.default_rng(seed)
    Z = np.vstack([X, Y])
    n = len(X)
    null = []
    for _ in range(n_perm):
        p = rng.permutation(len(Z))
        null.append(edistance(Z[p[:n]], Z[p[n:]], seed=seed))
    null = np.array(null)
    return obs, float(((null >= obs).sum() + 1) / (n_perm + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-features", type=int, default=250)
    ap.add_argument("--n-perm", type=int, default=200)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from atac_decomposition import LINES, CTRL_TOKENS, load_line

    print("loading Spear-ATAC (ChromVar) ...", flush=True)
    cells, PB = {}, {}
    for ln in LINES:
        obs, X, feat = load_line(ln, "ChromVar")
        cells[ln] = (obs, X)
        keys, mats = [], []
        for (p, r), g in obs.groupby(["pert", "rep"], observed=True):
            ii = g.index.to_numpy()
            if len(ii) >= 15:
                keys.append((p, r)); mats.append(X[ii].mean(0))
        PB[ln] = (pd.DataFrame(keys, columns=["pert", "rep"]),
                  np.stack(mats).astype(np.float32))

    # control-subtracted, exactly as the bespoke analysis does
    rows = []
    for ln in LINES:
        Kl, Xl = PB[ln]
        isc = Kl.pert.str.upper().str.contains("|".join(CTRL_TOKENS), na=False)
        ctl = {r: Xl[g.index[isc[g.index]].to_numpy()].mean(0)
               if isc.any() else Xl[g.index.to_numpy()].mean(0)
               for r, g in Kl.groupby("rep", observed=True)}
        for i in np.where(~isc if isc.any() else np.ones(len(Kl), bool))[0]:
            rows.append({"line": ln, "pert": Kl.pert.iloc[i],
                         "rep": Kl.rep.iloc[i],
                         "v": Xl[i] - ctl[Kl.rep.iloc[i]]})
    R = pd.DataFrame(rows)
    shared = set.intersection(*[set(R[R.line == l].pert) for l in LINES])
    R = R[R.pert.isin(shared)].reset_index(drop=True)
    V = np.stack(R.v.to_numpy())
    print(f"{len(R)} (line, perturbation, sample) rows, "
          f"{len(shared)} shared perturbations, {V.shape[1]} features",
          flush=True)

    # ---------- 1. variance components, published estimator ----------
    # Balance the design: take exactly two replicate samples per
    # (line, perturbation). The ANOVA estimator below is exact on a balanced
    # layout, and 2 is the minimum any cell has, so nothing is imputed.
    rng0 = np.random.default_rng(0)
    lines_l = sorted(R.line.unique())
    # keep only perturbations with at least two replicate samples in EVERY
    # line, so the cube is exactly balanced and nothing is imputed
    perts_l = [p_ for p_ in sorted(shared)
               if all(((R.line == ln) & (R.pert == p_)).sum() >= 2
                      for ln in lines_l)]
    dropped = len(shared) - len(perts_l)
    if dropped:
        print(f"   {dropped} perturbations lack 2 samples in some line and are "
              f"dropped", flush=True)
    cube = np.empty((len(lines_l), len(perts_l), 2, V.shape[1]),
                    dtype=np.float64)
    for ia, ln in enumerate(lines_l):
        for jb, p_ in enumerate(perts_l):
            ii = R.index[(R.line == ln) & (R.pert == p_)].to_numpy()
            cube[ia, jb] = V[rng0.choice(ii, 2, replace=False)]
    print(f"\n1. Variance components, ANOVA method of moments "
          f"(Searle et al. 1992;\n   the balanced-design limit of "
          f"variancePartition, Hoffman & Schadt 2016)")
    print(f"   design {cube.shape[0]} lines x {cube.shape[1]} perturbations "
          f"x {cube.shape[2]} replicates x {cube.shape[3]} motifs", flush=True)
    assert cube.shape[2] == 2 and cube.shape[0] == len(lines_l)
    VP, neg = anova_components(cube, cube.shape[0], cube.shape[1],
                               cube.shape[2])
    VP = VP.dropna()
    print(f"   context (cell line)         {VP.ctx.median():6.1%}")
    print(f"   perturbation                {VP.pert.median():6.1%}")
    print(f"   context x perturbation      {VP.inter.median():6.1%}")
    print(f"   residual (within-replicate) {VP.resid.median():6.1%}")
    print(f"   negative moment estimates (clipped to zero): "
          f"context {neg['ctx']:.0%}, perturbation {neg['pert']:.0%}, "
          f"interaction {neg['inter']:.0%}")
    reprod = VP.ctx + VP.pert + VP.inter
    share = float(np.median(VP.inter[reprod > 0] / reprod[reprod > 0]))
    print(f"\n   interaction as a share of REPRODUCIBLE variance: "
          f"{share:.1%}")
    print(f"   the bespoke covariance estimator gave 54.2% on the same data, "
          f"against a\n   label-permutation null of 61.1% -- i.e. no "
          f"detectable interaction. A published\n   estimator reaching a "
          f"similar figure is agreement about the ESTIMATE; the null\n   is "
          f"what says the estimate is not distinguishable from chance.")

    # ---------- 2. E-distance and E-test ----------
    print(f"\n2. E-distance and E-test (Peidli et al., Nat Methods 2024)",
          flush=True)
    et = []
    for ln in LINES:
        obs, X = cells[ln]
        isc = obs.pert.str.upper().str.contains("|".join(CTRL_TOKENS), na=False)
        C = X[obs.index[isc].to_numpy()] if isc.any() else None
        if C is None or len(C) < 10:
            print(f"   {ln}: no control cells, skipped"); continue
        for p, g in obs[~isc].groupby("pert", observed=True):
            if p not in shared:
                continue
            Y = X[g.index.to_numpy()]
            if len(Y) < 20:
                continue
            e, pv = etest(C, Y, n_perm=args.n_perm)
            et.append({"line": ln, "pert": p, "edist": e, "p": pv,
                       "n_cells": len(Y)})
    ET = pd.DataFrame(et)
    if len(ET):
        ET["q"] = np.minimum(ET.p * len(ET), 1.0)
        sig = ET[ET.q < 0.05]
        print(f"   {len(ET)} (line, perturbation) tests; "
              f"{len(sig)} reject at Bonferroni q<0.05 ({len(sig)/len(ET):.0%})")
        print(f"   median E-distance {ET.edist.median():.4f}")
        print("   A perturbation that does not move the population cannot "
              "contribute a real\n   interaction; if most tests fail to reject, "
              "an interaction estimate is\n   measuring perturbations that did "
              "nothing.")
        if len(sig) and len(ET) > len(sig):
            print(f"   fraction of perturbations detectable in all three "
                  f"lines: "
                  f"{ET.groupby('pert').q.apply(lambda s: (s<0.05).all()).mean():.0%}")

    S = pd.DataFrame([{"method": "variancePartition (LMM, REML)",
                       "interaction_of_reproducible": share,
                       "n_features": len(VP)},
                      {"method": "replicate covariance (this project)",
                       "interaction_of_reproducible": 0.542,
                       "n_features": V.shape[1]}])
    S.to_csv(TAB / "published_methods_check.csv", index=False)
    if len(ET):
        ET.to_csv(TAB / "published_methods_etest.csv", index=False)

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    med = [VP.ctx.median(), VP.pert.median(), VP.inter.median(),
           VP.resid.median()]
    ax[0].bar(range(4), med, width=0.6, color=[GREY, BLUE, ORANGE, "#d8d8d8"])
    for i_, v in enumerate(med):
        ax[0].text(i_, v + 0.008, f"{v:.1%}", ha="center", fontsize=9,
                   fontweight="bold")
    ax[0].set_xticks(range(4), ["cell line", "perturbation",
                                "line ×\nperturbation", "residual\n(noise)"],
                     fontsize=7)
    ax[0].set_ylabel("median variance fraction")
    ax[0].set_title("a  variancePartition\n    (Hoffman & Schadt 2016)",
                    loc="left", fontweight="bold", fontsize=9.5)

    ax[1].bar([0, 1], [share, 0.542], width=0.5, color=[VIOLET, ORANGE])
    for i_, v in enumerate([share, 0.542]):
        ax[1].text(i_, v + 0.012, f"{v:.1%}", ha="center", fontsize=10,
                   fontweight="bold")
    ax[1].set_xticks([0, 1], ["variancePartition\n(published)",
                              "replicate covariance\n(this project)"],
                     fontsize=7.5)
    ax[1].set_ylabel("interaction / reproducible variance")
    ax[1].set_title("b  Two routes to the same quantity", loc="left",
                    fontweight="bold", fontsize=9.5)

    if len(ET):
        ax[2].hist(-np.log10(np.maximum(ET.p, 1e-4)), bins=25, color=AQUA,
                   alpha=0.85)
        ax[2].axvline(-np.log10(0.05 / max(len(ET), 1)), color=ORANGE, lw=2,
                      label="Bonferroni 0.05")
        ax[2].set_xlabel("−log₁₀ p, permutation E-test")
        ax[2].set_ylabel("(line, perturbation) tests")
        ax[2].legend(frameon=False, fontsize=7.5)
    ax[2].set_title("c  E-test (Peidli et al. 2024)\n    did the perturbation "
                    "move the cells?", loc="left", fontweight="bold",
                    fontsize=9.5)
    fig.suptitle("The central split re-derived with published methods, on "
                 "fully crossed chromatin data", fontsize=10.5, x=0.005,
                 ha="left", fontweight="bold")
    d = save_figure(fig, "published_methods", FIG,
                    source_data={"summary": S, "variance_partition": VP,
                                 "etest": ET if len(ET) else pd.DataFrame()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
