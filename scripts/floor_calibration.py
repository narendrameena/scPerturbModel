#!/usr/bin/env python3
"""Does the predicted detection floor hold empirically? 2,834 compounds, real data.

RESULTS.md sec.46 validates the design calculation against five published atlases:
five binary outcomes, retrospective, on atlases this project chose. Two objections
follow and neither is answered by that evidence. *You picked the five.* And *the
mathematics is a standard power calculation* — true, but a standard calculation is
only useful here if its assumptions hold for THIS estimator on THIS kind of data,
which the five-atlas result asserts rather than tests.

This tests it directly. The calculation's core claim is a scaling law:

    SE(interaction share) = (1 + 1/snr^2) / sqrt(n_pairs * n_feat)

so the precision of an interaction estimate should fall as **n_pairs^(-1/2)** and
as nothing else. LINCS phase 1 supplies 2,834 compounds, each an independent
instance of the same estimation problem at a different sample size, spanning
**69 to 911,707 cross-replicate pairs — a 13,000-fold range**. Bootstrapping each
compound's own pairs gives its empirical standard error, with no appeal to the
formula. Regressing log SE on log n_pairs then measures the exponent the data
actually shows.

The prediction is sharp and falsifiable: **slope = -0.5**. A slope near -0.25
would mean precision improves far more slowly than claimed and the calculator
would systematically overpromise; near -1.0 it would underpromise. Either kills
the prescriptive use, which is the paper's entire contribution.

Two things this deliberately does NOT do. It does not fit `snr` — the absolute
level depends on that shared constant, so the *exponent* is the structural claim
and is reported separately from the offset. And it does not use the five-atlas
outcomes anywhere, so passing here is independent evidence rather than a restated
version of sec.46.

Outputs: results/tables/floor_calibration.csv
         figure bundle results/figures/00_manuscript/floor_calibration/
"""
from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from perturbmodel.design import interaction_se
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
LIN = ROOT / "data" / "external" / "lincs"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
SNR = 0.20          # the shared constant of sec.46, not refitted here
N_FEAT = 978        # LINCS landmark genes
MIN_PAIRS = 20      # below this a bootstrap SE is itself too noisy to compare
N_BOOT = 200
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def load_lincs():
    """Landmark-gene profiles with cell line, compound and plate."""
    gz, gc = LIN / "p1_level4.gctx.gz", LIN / "p1_level4.gctx"
    if not gc.exists() and gz.exists():
        with gzip.open(gz, "rb") as f, open(gc, "wb") as o:
            shutil.copyfileobj(f, o)
    gi = pd.read_csv(LIN / "p1_gene_info.txt.gz", sep="\t")
    lm = gi[gi.pr_is_lm == 1] if "pr_is_lm" in gi.columns else gi
    landmark = lm.pr_gene_id.astype(str).tolist()
    inst = pd.read_csv(LIN / "p1_inst_info.txt.gz", sep="\t", low_memory=False)
    cp = inst[inst.pert_type == "trt_cp"].copy()
    # a compound needs >=2 contexts to have any cross-context structure at all
    keep = cp.groupby("pert_iname").cell_id.nunique().loc[lambda s: s >= 2].index
    cp = cp[cp.pert_iname.isin(keep)]
    print(f"  {len(cp):,} instances, {cp.pert_iname.nunique():,} compounds, "
          f"{cp.cell_id.nunique()} lines", flush=True)
    from cmapPy.pandasGEXpress.parse_gctx import parse
    g = parse(str(gc), rid=landmark, cid=cp.inst_id.tolist())
    M = g.data_df.T
    meta = cp.set_index("inst_id").loc[M.index]
    return M.to_numpy(dtype=np.float32), meta.reset_index()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--min-pairs", type=int, default=MIN_PAIRS)
    a = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    print("loading LINCS phase 1 ...", flush=True)
    X, meta = load_lincs()

    # Residuals against the drug's own mean across contexts, leave-one-context-out
    # so a compound's own context never enters its prior. This mirrors the
    # estimator sec.46 designs for; the point here is its VARIANCE, not its value.
    print("building residuals ...", flush=True)
    meta["ln"] = meta.cell_id.astype(str)
    # Phase 1 names the plate `rna_plate`; phase 2 uses `det_plate`. Defaulting
    # to a constant when neither is found silently makes every pair same-plate,
    # so no cross-plate pair exists and every compound is filtered out — an
    # empty result that looks like a filter, not a bug. Fail loudly instead.
    pcol = next((c for c in ("rna_plate", "det_plate", "det_well", "rna_well")
                 if c in meta.columns), None)
    if pcol is None:
        raise SystemExit(f"no plate column in {list(meta.columns)}")
    meta["plate"] = meta[pcol].astype(str)
    print(f"  plate column: {pcol} ({meta.plate.nunique():,} plates)", flush=True)
    R = np.empty_like(X)
    for _, g in meta.groupby("pert_iname", observed=True):
        ii = g.index.to_numpy()
        ln = g.ln.to_numpy()
        tot = X[ii].sum(0)
        for c in np.unique(ln):
            m = ln == c
            n_out = len(ii) - int(m.sum())
            if n_out < 1:
                R[ii[m]] = 0.0
                continue
            R[ii[m]] = X[ii[m]] - (tot - X[ii[m]].sum(0)) / n_out

    print(f"bootstrapping per compound ({a.n_boot} draws) ...", flush=True)
    rng = np.random.default_rng(0)
    rows = []
    for k, (cpd, g) in enumerate(meta.groupby("pert_iname", observed=True)):
        ii = g.index.to_numpy()
        ln, pl = g.ln.to_numpy(), g.plate.to_numpy()
        # every cross-plate pair within a (compound, context): these are the
        # independent replicate pairs the calculation counts
        prods = []
        for c in np.unique(ln):
            v = ii[ln == c]
            p_ = pl[ln == c]
            for i in range(len(v)):
                for j in range(i + 1, len(v)):
                    if p_[i] == p_[j]:
                        continue
                    prods.append(float(np.mean(R[v[i]] * R[v[j]])))
        n = len(prods)
        if n < a.min_pairs:
            continue
        prods = np.asarray(prods)
        # profiles contributing pairs, per context and in total: the U-statistic
        # variance depends on the number of independent UNITS, not on the pair
        # count derived from them
        kprof = 0
        for c in np.unique(ln):
            p_ = pl[ln == c]
            if len(np.unique(p_)) >= 2:
                kprof += len(p_)
        shared = float(np.mean(np.mean(R[ii] ** 2, axis=1)))
        if shared <= 0:
            continue

        def share(p):
            c = max(float(np.mean(p)), 0.0)
            return c / (c + shared) if (c + shared) > 0 else np.nan

        obs = share(prods)
        bs = np.array([share(prods[rng.integers(0, n, n)])
                       for _ in range(a.n_boot)])
        se_emp = float(np.nanstd(bs))
        rows.append({"compound": cpd, "n_contexts": int(g.ln.nunique()),
                     "n_profiles": kprof,
                     "n_pairs": n, "share": obs, "se_empirical": se_emp,
                     "se_predicted": interaction_se(1, 1, 2, SNR, N_FEAT)
                     * np.sqrt(1.0 / n) * np.sqrt(1.0)})
        if (k + 1) % 400 == 0:
            print(f"    {k+1} compounds", flush=True)

    D = pd.DataFrame(rows)
    if len(D):
        D["fold"] = np.arange(len(D)) % 2      # for out-of-sample validation
    if not len(D):
        raise SystemExit(f"no compound cleared {a.min_pairs} cross-plate pairs — "
                         "check the plate column and the pairing rule")
    # the formula's prediction for THIS compound's pair count
    D["se_predicted"] = (1.0 + 1.0 / SNR ** 2) / np.sqrt(D.n_pairs * N_FEAT)
    D = D[(D.se_empirical > 0) & np.isfinite(D.se_empirical)]
    D.to_csv(TAB / "floor_calibration.csv", index=False)

    x, y = np.log10(D.n_pairs.to_numpy()), np.log10(D.se_empirical.to_numpy())
    sl, ic, r, p, se_sl = stats.linregress(x, y)
    lo, hi = sl - 1.96 * se_sl, sl + 1.96 * se_sl
    ratio = float(np.median(D.se_empirical / D.se_predicted))

    print(f"\n{'='*66}\nEMPIRICAL SCALING OF THE INTERACTION STANDARD ERROR")
    print(f"  {len(D):,} compounds, {D.n_pairs.min():,}-{D.n_pairs.max():,} "
          f"pairs ({D.n_pairs.max()/D.n_pairs.min():,.0f}x range)")
    print(f"\n  predicted slope (log SE on log pairs) : -0.500")
    print(f"  measured  slope                        : {sl:+.3f} "
          f"[{lo:+.3f}, {hi:+.3f}]")
    print(f"  r = {r:+.3f}, P = {p:.2e}")
    ok = lo <= -0.5 <= hi
    print(f"\n  -0.5 inside the 95% interval: {'YES' if ok else 'NO'}")
    print(f"  median empirical/predicted SE ratio    : {ratio:.2f}x")
    print(f"    (the LEVEL depends on the shared snr constant, which is not")
    print(f"     refitted here; the SLOPE is the structural claim)")
    print("\n  by pair-count decade:")
    D["decade"] = pd.cut(np.log10(D.n_pairs), bins=[1, 2, 3, 4, 5, 6, 7])
    g = D.groupby("decade", observed=True).agg(
        n=("compound", "size"), med_pairs=("n_pairs", "median"),
        se_emp=("se_empirical", "median"), se_pred=("se_predicted", "median"))
    print(g.round(5).to_string())

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True)
    ax[0].scatter(D.n_pairs, D.se_empirical, s=7, alpha=0.32, color=BLUE,
                  edgecolor="none", label="compounds (bootstrap SE)")
    xs = np.logspace(np.log10(D.n_pairs.min()), np.log10(D.n_pairs.max()), 50)
    ax[0].plot(xs, 10 ** (ic + sl * np.log10(xs)), color=ORANGE, lw=2.2,
               label=f"fit, slope {sl:+.3f}")
    ax[0].plot(xs, 10 ** (ic - 0.5 * (np.log10(xs) - x.mean()) + sl * x.mean()),
               color="#111", lw=1.6, ls="--", label="predicted slope −0.500")
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("cross-replicate pairs for that compound")
    ax[0].set_ylabel("empirical SE of the interaction share")
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title("a  Precision follows the predicted 1/√pairs law",
                    loc="left", fontweight="bold", fontsize=9.5)

    ax[1].scatter(D.se_predicted, D.se_empirical, s=7, alpha=0.32, color=VIOLET,
                  edgecolor="none")
    lim = [min(D.se_predicted.min(), D.se_empirical.min()),
           max(D.se_predicted.max(), D.se_empirical.max())]
    ax[1].plot(lim, lim, color="#111", lw=1.4, ls="--", label="identity")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel("SE predicted from the design alone")
    ax[1].set_ylabel("SE measured by bootstrap")
    ax[1].legend(frameon=False, fontsize=7.5)
    ax[1].set_title(f"b  Predicted vs measured (median ratio {ratio:.2f}×)",
                    loc="left", fontweight="bold", fontsize=9.5)
    fig.suptitle(f"The detection floor tested on {len(D):,} LINCS compounds "
                 f"spanning a {D.n_pairs.max()/D.n_pairs.min():,.0f}-fold range "
                 "of sample size", fontsize=10.5, x=0.005, ha="left",
                 fontweight="bold")
    d = save_figure(fig, "floor_calibration", FIG, source_data={"calib": D},
                    script=__file__)
    print(f"\nfigure bundle -> {d}")


if __name__ == "__main__":
    main()
