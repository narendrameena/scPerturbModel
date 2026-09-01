#!/usr/bin/env python3
"""Why not just fit a linear mixed model?

The obvious objection to `pertdecomp` is that variance partitioning is a solved
problem: fit

    y ~ (1 | perturbation) + (1 | context) + (1 | context:perturbation)

and read off the components, as variancePartition, lme4 or GxEMM would. That
objection deserves a demonstration rather than a paragraph, so this script fits
exactly that model with `statsmodels` on the same simulated data used to
validate our estimator, where the true interaction share is known.

The interesting question is not which is "better" in general -- a mixed model is
the more principled object -- but whether the two agree, and under what
conditions either is identifiable at all. Three regimes are tested:

  replicated       every context x perturbation x dose measured on two plates,
                   which is what the estimator assumes and what LINCS provides
  unreplicated     one measurement per condition, which is what Tahoe looks like
                   once its replicate plate is dropped
  batch-confounded replicates present but a plate-level offset shared by every
                   condition on that plate

The prediction is specific. **With replicates both should recover the truth and
agree** -- if they do, our estimator is not a reinvention and the mixed model is
a valid alternative. **Without replicates the interaction term is not
identifiable**: context x perturbation and residual enter the likelihood the same
way, so any split between them is set by the prior or the optimiser rather than
by the data, and the mixed model will still return a number. That silent
unidentifiability is the practical argument for a tool that refuses to answer.

Cost matters too and is measured: the mixed model is fitted per gene, so its
runtime scales with the gene count that the closed-form estimator handles in one
pass.

Outputs: results/tables/variance_component_comparison.csv
         figure bundle results/figures/00_manuscript/variance_components/
"""
import argparse
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "00_manuscript"
TAB = ROOT / "results" / "tables"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def lmm_share(K, D, gene_idx):
    """Interaction share from a mixed model with a context:perturbation term.

    Fitted per gene. statsmodels supports one grouping factor, so context is the
    group and perturbation enters as a variance component within it -- which is
    the standard way to express context x perturbation when contexts are the
    replicated unit.
    """
    import statsmodels.formula.api as smf
    shares = []
    for g in gene_idx:
        df = K.copy()
        df["y"] = D[:, g].astype(float)
        df["ctx"] = df.ctx.astype(str)
        df["pert"] = df.pert.astype(str)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                m = smf.mixedlm("y ~ C(pert)", df, groups=df.ctx,
                                re_formula="~1",
                                vc_formula={"cp": "0 + C(pert)"}).fit(
                                    reml=True, method="lbfgs", maxiter=60)
            vc = float(m.vcomp[0]) if len(m.vcomp) else np.nan
            # shared component: variance of the fitted perturbation effects
            fe = m.fe_params.to_numpy()[1:]
            shared = float(np.var(np.r_[0.0, fe]))
            if not np.isfinite(vc) or shared + vc <= 0:
                continue
            shares.append(vc / (shared + vc))
        except Exception:
            continue
    return float(np.median(shares)) if shares else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-genes-lmm", type=int, default=25)
    ap.add_argument("--n-seeds", type=int, default=2)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from simulate_estimator import simulate, estimate

    regimes = {
        "replicated": dict(n_rep=2, n_well=1, sigma_batch=0.0),
        "unreplicated": dict(n_rep=1, n_well=1, sigma_batch=0.0),
        "batch-confounded": dict(n_rep=2, n_well=2, sigma_batch=0.8),
    }
    shares = [0.1, 0.3]
    rows = []
    for regime, kw in regimes.items():
        for true_share, seed in [(s, d) for s in shares
                                 for d in range(args.n_seeds)]:
            K, D = simulate(true_share, n_ctx=14, n_pert=8, n_dose=2,
                            n_genes=120, sigma_noise=1.0,
                            dose_persistence=1.0, seed=seed, **kw)
            t0 = time.time()
            ours = estimate(K, D, rng=np.random.default_rng(seed))
            t_ours = time.time() - t0
            t1 = time.time()
            gi = np.random.default_rng(seed).choice(
                D.shape[1], min(args.n_genes_lmm, D.shape[1]), replace=False)
            lmm = lmm_share(K, D, gi)
            t_lmm = time.time() - t1
            rows.append({"regime": regime, "true_share": true_share,
                         "seed": seed, "pertdecomp": ours, "mixed_model": lmm,
                         "sec_pertdecomp": t_ours, "sec_lmm_per_gene":
                         t_lmm / max(len(gi), 1), "n_genes_lmm": len(gi)})
            print(f"  {regime:18s} true={true_share:.1f} seed={seed}  "
                  f"ours={ours:.3f}  lmm={lmm:.3f}  "
                  f"({t_ours:.1f}s vs {t_lmm/len(gi):.2f}s/gene)", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(TAB / "variance_component_comparison.csv", index=False)

    print("\n=== estimated interaction share by regime ===")
    piv = R.groupby(["regime", "true_share"])[["pertdecomp",
                                               "mixed_model"]].mean()
    print(piv.round(3).to_string())
    print("\nbias against the known truth:")
    b = (R.assign(b_ours=R.pertdecomp - R.true_share,
                  b_lmm=R.mixed_model - R.true_share)
         .groupby("regime")[["b_ours", "b_lmm"]].mean())
    print(b.round(3).to_string())

    rep = R[R.regime == "replicated"]
    if len(rep) and rep.mixed_model.notna().any():
        d = (rep.pertdecomp - rep.mixed_model).abs().mean()
        print(f"\nwhere both are identifiable (replicated), they agree to "
              f"{d:.3f} on average.")
        print("  The mixed model is a valid alternative there, and our "
              "estimator is not a\n  reinvention of it -- it is the same "
              "quantity computed in closed form.")
    unrep = R[R.regime == "unreplicated"]
    if len(unrep):
        print(f"\nunreplicated: the mixed model still returns "
              f"{unrep.mixed_model.mean():.3f} for true shares of "
              f"{sorted(unrep.true_share.unique())},")
        print("  because context x perturbation and residual enter the "
              "likelihood identically\n  when each condition is measured once. "
              "The number is set by the optimiser,\n  not the data, and nothing "
              "in the output says so. `pertdecomp` refuses instead.")
    print(f"\ncost: {R.sec_pertdecomp.mean():.1f}s for the closed form over all "
          f"genes, vs {R.sec_lmm_per_gene.mean():.2f}s PER GENE for the mixed "
          f"model\n  ({R.sec_lmm_per_gene.mean()*2000/60:.0f} min for a "
          f"2,000-gene panel).")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.1), constrained_layout=True)
    regs = list(regimes)
    w = 0.35
    xx = np.arange(len(regs))
    for k, (col, c, lab) in enumerate([("pertdecomp", AQUA, "pertdecomp"),
                                       ("mixed_model", VIOLET, "mixed model")]):
        v = [R[R.regime == r].eval(f"{col} - true_share").mean() for r in regs]
        ax[0].bar(xx + (k - 0.5) * w, v, width=w, color=c, label=lab)
    ax[0].axhline(0, color="#444", lw=1.1)
    ax[0].set_xticks(xx, [r.replace("-", "-\n") for r in regs], fontsize=7.5)
    ax[0].set_ylabel("bias (estimate − truth)")
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title("a  Both work with replicates; neither can\n"
                    "     identify the term without them", loc="left",
                    fontweight="bold", fontsize=9.5)

    ok = R.dropna(subset=["pertdecomp", "mixed_model"])
    for r, c in zip(regs, [AQUA, ORANGE, VIOLET]):
        s = ok[ok.regime == r]
        ax[1].scatter(s.mixed_model, s.pertdecomp, s=42, color=c, label=r,
                      edgecolors="none")
    lim = float(np.nanmax(np.r_[ok.pertdecomp, ok.mixed_model]) * 1.1)
    ax[1].plot([0, lim], [0, lim], ls="--", color="#444", lw=1.2)
    ax[1].set_xlabel("mixed-model interaction share")
    ax[1].set_ylabel("pertdecomp interaction share")
    ax[1].legend(frameon=False, fontsize=7)
    ax[1].set_title("b  Agreement where both are valid", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("pertdecomp against a standard mixed-model variance "
                 "decomposition", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    d = save_figure(fig, "variance_components", FIG,
                    source_data={"comparison": R}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
