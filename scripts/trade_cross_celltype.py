#!/usr/bin/env python3
"""The test §28 left standing: does the inflation appear in a REPLICATE-VALIDATED
context-specificity statistic published by another group?

§27 found our interaction shares inflated ~1.5x by a cell property left inside
the residual. §28 tried to show the same in other groups' published numbers and
failed: neither sci-Plex's 48% nor CMap's 26% reproduced under our
implementation, and where our versions did move they moved by 1.01-1.15x. The
diagnosis was that those statistics are computed on raw, noise-dominated
residuals, where a term worth a few percent of total variance changes nothing.
Our estimator instead compares components that noise cannot inflate -- each is a
covariance between two independent estimates -- and there the cell property is
about a third of the non-drug signal.

That left a prediction with no data to test it: **the inflation should appear in
noise-corrected statistics and not in raw ones.** Nadig, Replogle, Pogson et al.
(*Nature Genetics* 57:1228, 2025) provide exactly the missing case. TRADE splits
perturbation effects into a cross-cell-type-consistent and a cell-type-dependent
part -- 56%/44% -- and, decisively for us, their deposit ships **lfcSE beside
every log2 fold change** for four cell lines given the same essential-gene
perturbations (K562, RPE1, HepG2, Jurkat). The standard errors make the
noise-corrected version constructible, and having four cell lines makes the
consistent component estimable as a between-line covariance, which no amount of
measurement noise can inflate because two cell lines' noise is independent.

Both statistics are computed on identical data, and the same cell-line main
effect is removed from both:

  raw            consistent = mean(LFC across lines)^2
                 total      = mean(LFC^2 across lines)
                 -- noise enters both terms, and enters `total` twice as hard

  noise-corrected consistent = <LFC in line a, LFC in line b> over line pairs
                 reproducible = mean over lines of var(LFC) - mean(lfcSE^2)
                 -- the first is noise-free by independence, the second has the
                    sampling variance subtracted using the published SEs

The correction removed is each line's mean LFC across the OTHER perturbations,
leave-one-perturbation-out, so a genuine perturbation-specific effect cannot be
absorbed into it.

If the raw split barely moves and the noise-corrected split moves substantially,
the §28 prediction holds and the defect is a property of noise-corrected
estimators generally -- which is a statement about the method, not about us.

Outputs: results/tables/trade_cross_celltype.csv
         figure bundle results/figures/00_manuscript/trade_cross_celltype/
"""
import argparse
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
TD = ROOT / "data" / "external" / "trade"
FIG = ROOT / "results" / "figures" / "00_manuscript"
TAB = ROOT / "results" / "tables"
LINES = ("K562", "RPE1", "HepG2", "Jurkat")
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"


def load(line):
    """LFC and lfcSE for one cell line, genes x perturbations."""
    L = pd.read_csv(TD / f"{line}Essential_Log2FoldChange.csv.gz", index_col=0)
    S = pd.read_csv(TD / f"{line}Essential_lfcSE.csv.gz", index_col=0)
    g = L.index.intersection(S.index)
    p = L.columns.intersection(S.columns)
    print(f"  {line:7s} {L.shape} -> {len(g)} genes x {len(p)} perturbations",
          flush=True)
    return L.loc[g, p], S.loc[g, p]


def loo_line_effect(M, min_perts=20):
    """Each line's mean LFC across the OTHER perturbations, CENTRED over lines.

    Two properties are required and the second is easy to omit:

      * leave-one-perturbation-out, so the correction cannot absorb the
        perturbation-specific effect it is meant to leave alone;
      * **centred across cell lines**, so the correction is a CONTRAST between
        lines rather than an absolute mean.

    Without the centring the correction is not a line property at all. Every
    perturbation in this deposit is an essential-gene knockdown, and knocking
    down any essential gene drives a common growth-arrest programme; each line's
    uncentred mean is dominated by that shared programme, so subtracting it
    removes signal the statistic is supposed to keep. Measured on this data the
    four uncentred vectors correlate with each other at r = +0.51 and with the
    cross-line shared response at r = +0.80 -- they are one programme, not four
    line properties. ``perturbmodel.celldrug`` centres for exactly this reason;
    this function did not, and the first version of this analysis was measuring
    the omission rather than the correction.

    Takes the full lines x genes x perturbations tensor so the centring is
    possible, and returns a correction of the same shape.
    """
    n = M.shape[2]
    if n < min_perts:
        return np.zeros_like(M)
    tot = np.nansum(M, axis=2, keepdims=True)
    cnt = np.sum(np.isfinite(M), axis=2, keepdims=True)
    A = (tot - np.nan_to_num(M)) / np.maximum(cnt - np.isfinite(M), 1)
    return A - np.nanmean(A, axis=0, keepdims=True)


def split_raw(M):
    """The split as a raw variance ratio, noise included in both terms."""
    mu = np.nanmean(M, axis=0)                    # over lines
    consistent = mu ** 2
    total = np.nanmean(M ** 2, axis=0)
    ok = np.isfinite(consistent) & np.isfinite(total) & (total > 0)
    return (float(np.nansum(consistent[ok]) / np.nansum(total[ok])),
            consistent, total, ok)


def split_corrected(M, SE):
    """The split with noise removed from both terms.

    consistent   between-line covariance, averaged over the six line pairs.
                 Two cell lines are measured in separate experiments, so their
                 sampling noise is independent and contributes zero in
                 expectation -- the same device the project's own estimator uses
                 with replicate plates.
    reproducible per-line variance minus the mean published sampling variance,
                 so what remains is the part of the effect that would survive
                 remeasurement.
    """
    nL = M.shape[0]
    covs = []
    for a, b in combinations(range(nL), 2):
        covs.append(M[a] * M[b])
    consistent = np.nanmean(np.stack(covs), axis=0)
    per_line_var = np.nanmean(M ** 2, axis=0)
    noise = np.nanmean(SE ** 2, axis=0)
    reproducible = per_line_var - noise
    ok = (np.isfinite(consistent) & np.isfinite(reproducible)
          & (reproducible > 0))
    return (float(np.nansum(consistent[ok]) / np.nansum(reproducible[ok])),
            consistent, reproducible, ok, float(np.nansum(noise[ok])
                                                / np.nansum(per_line_var[ok])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--published-consistent", type=float, default=0.56)
    ap.add_argument("--max-perts", type=int, default=0,
                    help="0 = use every shared perturbation")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    print("loading TRADE deposit (Nadig et al. 2025) ...", flush=True)
    D = {}
    for ln in LINES:
        D[ln] = load(ln)
    genes = D[LINES[0]][0].index
    perts = D[LINES[0]][0].columns
    for ln in LINES[1:]:
        genes = genes.intersection(D[ln][0].index)
        perts = perts.intersection(D[ln][0].columns)
    perts = list(perts)
    if args.max_perts:
        perts = perts[:args.max_perts]
    print(f"\nshared across all {len(LINES)} lines: {len(genes)} genes x "
          f"{len(perts)} perturbations", flush=True)
    if len(perts) < 50 or len(genes) < 500:
        print("too little overlap to test"); return

    M = np.stack([D[ln][0].loc[genes, perts].to_numpy(np.float32)
                  for ln in LINES])          # lines x genes x perts
    SE = np.stack([D[ln][1].loc[genes, perts].to_numpy(np.float32)
                   for ln in LINES])
    del D
    print(f"tensor {M.shape}", flush=True)

    # ---- is the "line main effect" actually a line property here? --------
    # The correction assumes a context's mean response across other
    # perturbations is a property of that context. That holds only if the
    # perturbations push in diverse directions. Every perturbation in this
    # deposit is an ESSENTIAL-gene knockdown, and knocking down any essential
    # gene produces a common growth-arrest signature -- so the per-line mean may
    # be dominated by that SHARED programme rather than by anything
    # line-specific. Two checks settle it:
    #   * if the four lines' correction vectors correlate with each other, they
    #     are capturing something shared, not four separate line properties;
    #   * if each correlates with the cross-line mean response, the correction
    #     is removing signal the statistic is supposed to keep.
    ALPHA = np.stack([np.nanmean(M[i], axis=1) for i in range(len(LINES))])
    ALPHA_C = ALPHA - np.nanmean(ALPHA, axis=0, keepdims=True)
    shared = np.nanmean(np.nanmean(M, axis=0), axis=1)
    ok_a = np.isfinite(ALPHA).all(0) & np.isfinite(shared)
    pair_r = []
    for a, b in combinations(range(len(LINES)), 2):
        pair_r.append(float(np.corrcoef(ALPHA[a][ok_a], ALPHA[b][ok_a])[0, 1]))
    with_shared = [float(np.corrcoef(ALPHA[i][ok_a], shared[ok_a])[0, 1])
                   for i in range(len(LINES))]
    pair_rc = [float(np.corrcoef(ALPHA_C[a][ok_a], ALPHA_C[b][ok_a])[0, 1])
               for a, b in combinations(range(len(LINES)), 2)]
    shared_rc = [float(np.corrcoef(ALPHA_C[i][ok_a], shared[ok_a])[0, 1])
                 for i in range(len(LINES))]
    print(f"\nIS THE CORRECTION A LINE PROPERTY HERE?")
    print(f"  correlation between different lines' correction vectors: "
          f"median r = {np.median(pair_r):+.3f}  {[round(x,2) for x in pair_r]}")
    print(f"  correlation of each with the cross-line shared response:  "
          f"median r = {np.median(with_shared):+.3f}  "
          f"{[round(x,2) for x in with_shared]}")
    print(f"  after centring across lines (what is actually "
          f"subtracted):")
    print(f"    between lines r = {np.median(pair_rc):+.3f}   "
          f"vs shared response r = {np.median(shared_rc):+.3f}")
    if np.median(pair_rc) > 0.5:
        print("  The four lines' 'line effects' are strongly correlated WITH "
              "EACH OTHER, so\n  they are not four line properties -- they are "
              "one shared programme. Every\n  perturbation here is an "
              "essential-gene knockdown, which drives a common\n  "
              "growth-arrest response, and the leave-one-perturbation-out mean "
              "absorbs it.\n  The correction is CONFOUNDED in this dataset and "
              "its result should not be\n  read as evidence either way.")
    else:
        print("  The correction vectors are largely independent between lines, "
              "so they do\n  behave as line properties and the comparison "
              "below is interpretable.")

    rows = []
    for tag, Mx in (("as published (line effect in)", M),
                    ("cell-line main effect removed",
                     M - loo_line_effect(M))):
        r_raw, c_r, t_r, ok_r = split_raw(Mx)
        r_cor, c_c, t_c, ok_c, noise_share = split_corrected(Mx, SE)
        print(f"\n{tag}")
        print(f"  RAW            consistent {r_raw:6.1%}  "
              f"cell-type-dependent {1-r_raw:6.1%}   "
              f"(n = {int(ok_r.sum()):,} gene x perturbation)")
        print(f"  NOISE-CORRECTED consistent {r_cor:6.1%}  "
              f"cell-type-dependent {1-r_cor:6.1%}   "
              f"(n = {int(ok_c.sum()):,})")
        print(f"  sampling noise is {noise_share:.1%} of the raw per-line "
              f"variance")
        rows.append({"variant": tag,
                     "alpha_between_line_r": float(np.median(pair_r)),
                     "alpha_vs_shared_r": float(np.median(with_shared)),
                     "consistent_raw": r_raw,
                     "dependent_raw": 1 - r_raw,
                     "consistent_corrected": r_cor,
                     "dependent_corrected": 1 - r_cor,
                     "noise_share_of_raw_variance": noise_share,
                     "n_raw": int(ok_r.sum()), "n_corrected": int(ok_c.sum())})

    T = pd.DataFrame(rows)
    T.to_csv(TAB / "trade_cross_celltype.csv", index=False)

    pub = args.published_consistent
    base = T.iloc[0]
    corr = T.iloc[1]
    print("\n" + "=" * 74)
    print(f"published (Nadig et al. 2025):          consistent {pub:.0%}  "
          f"cell-type-dependent {1-pub:.0%}")
    print(f"our replication, RAW:                   consistent "
          f"{base.consistent_raw:.0%}   -- "
          f"{'reproduced' if abs(base.consistent_raw - pub) < 0.10 else 'NOT reproduced'}")
    print(f"our replication, NOISE-CORRECTED:       consistent "
          f"{base.consistent_corrected:.0%}")
    d_raw = 100 * (base.dependent_raw - corr.dependent_raw)
    d_cor = 100 * (base.dependent_corrected - corr.dependent_corrected)
    print(f"\nremoving the cell-line main effect moves "
          f"'cell-type-dependent' (positive = the omission was inflating it)")
    print(f"  RAW statistic             {base.dependent_raw:6.1%} -> "
          f"{corr.dependent_raw:6.1%}   ({d_raw:+.1f} points)")
    print(f"  NOISE-CORRECTED statistic {base.dependent_corrected:6.1%} -> "
          f"{corr.dependent_corrected:6.1%}   ({d_cor:+.1f} points)")
    if d_raw < 0 and d_cor < 0:
        print("  Both are NEGATIVE — the correction makes the response look "
              "MORE\n  cell-type-dependent, the opposite of its effect in "
              "PRISM. Combined with the\n  confound check above, that is what "
              "a correction removing shared signal\n  rather than a context "
              "property looks like.")
    if abs(base.consistent_raw - pub) >= 0.10:
        print("\nThe raw replication does not match the published value, so no "
              "claim is made\nabout Nadig et al.'s number. The RAW-vs-CORRECTED "
              "contrast below is internal\nand still tests the §28 prediction, "
              "because both are computed here on the\nsame data with the same "
              "correction.")
    if abs(d_cor) > 2 * max(abs(d_raw), 1e-9) and np.median(pair_rc) <= 0.5:
        print("\nPREDICTION HOLDS: the cell-line main effect moves the "
              "noise-corrected split\nseveral times more than the raw one. The "
              "inflation is a property of\nnoise-corrected context-specificity "
              "estimators, which is a statement about\nthe method rather than "
              "about this project.")
    else:
        print("\nPREDICTION FAILS: the noise-corrected split does not move "
              "appreciably more\nthan the raw one. §28's explanation for why "
              "other groups' numbers were\nunaffected is then wrong too, and "
              "the defect should be reported as ours\nalone with no proposed "
              "generality.")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)

    xx = np.arange(2)
    w = 0.36
    ax[0].bar(xx - w / 2, [base.dependent_raw, base.dependent_corrected], w,
              color=ORANGE, label="line effect left in (as published)")
    ax[0].bar(xx + w / 2, [corr.dependent_raw, corr.dependent_corrected], w,
              color=AQUA, label="cell-line main effect removed")
    for i_, (a_, b_) in enumerate([(base.dependent_raw, corr.dependent_raw),
                                   (base.dependent_corrected,
                                    corr.dependent_corrected)]):
        ax[0].text(i_ - w / 2, a_ + .012, f"{a_:.0%}", ha="center", fontsize=8.5,
                   fontweight="bold")
        ax[0].text(i_ + w / 2, b_ + .012, f"{b_:.0%}", ha="center", fontsize=8.5,
                   fontweight="bold")
    ax[0].axhline(1 - pub, ls="--", color="#555", lw=1.2,
                  label=f"published {1-pub:.0%}")
    ax[0].set_xticks(xx, ["raw\n(noise included)",
                          "noise-corrected\n(covariance + published SEs)"],
                     fontsize=8)
    ax[0].set_ylabel("'cell-type-dependent' share")
    ax[0].legend(frameon=False, fontsize=7)
    ax[0].set_title("a  TRADE's split, recomputed two ways", loc="left",
                    fontweight="bold", fontsize=9.5)

    ax[1].bar([0, 1], [abs(d_raw), abs(d_cor)], width=0.5,
              color=[GREY, VIOLET])
    for i_, v in enumerate([abs(d_raw), abs(d_cor)]):
        ax[1].text(i_, v + 0.08, f"{v:.1f} pts", ha="center", fontsize=10,
                   fontweight="bold")
    ax[1].set_xticks([0, 1], ["raw", "noise-corrected"], fontsize=8.5)
    ax[1].set_ylabel("shift in the statistic (percentage points)")
    ax[1].set_title("b  The prediction from §28", loc="left",
                    fontweight="bold", fontsize=9.5)
    ax[1].text(0.5, 0.9, "§28 predicted the right bar\nis the larger one",
               transform=ax[1].transAxes, ha="center", fontsize=7,
               color="#555")

    ax[2].bar([0], [base.noise_share_of_raw_variance], width=0.4, color=ORANGE)
    ax[2].set_xticks([0], ["TRADE, 4 lines"], fontsize=8.5)
    ax[2].set_ylim(0, 1.0)
    ax[2].set_ylabel("sampling noise as a share of raw per-line variance")
    ax[2].text(0, base.noise_share_of_raw_variance + 0.02,
               f"{base.noise_share_of_raw_variance:.0%}", ha="center",
               fontsize=11, fontweight="bold")
    ax[2].text(0.5, 0.55, "This is why a raw statistic barely moves:\nmost of "
               "what it measures is noise, so removing\na real term worth a few "
               "percent changes little.", transform=ax[2].transAxes,
               ha="center", fontsize=7, color="#444")
    ax[2].set_title("c  Why raw statistics are insensitive", loc="left",
                    fontweight="bold", fontsize=9.5)
    fig.suptitle("Testing the §28 prediction on another group's "
                 "replicate-validated statistic (TRADE, Nadig et al. 2025)",
                 fontsize=10.5, x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "trade_cross_celltype", FIG, source_data={"summary": T},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
