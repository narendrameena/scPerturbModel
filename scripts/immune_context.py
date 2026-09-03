#!/usr/bin/env python3
"""Is drug response more or less context-dependent in immune contexts?

Everything so far has been dominated by solid-tumour lines, and the question of
whether immune contexts behave the same way matters for two reasons: most
therapeutic interest in context-specific response is now in immune cells, and
the OP3 result (RESULTS.md §9) hinted that immune contexts might differ but was
too underpowered per compound to say.

PRISM annotates lineage, giving 79 haematopoietic and lymphoid lines (52
leukaemia, 27 lymphoma) against ~660 solid lines, screened with the same
compounds on the same plates -- so the comparison is internal, with assay,
compound set and batch structure held fixed. That is a much cleaner contrast
than comparing two studies.

**The control that makes this interpretable.** We established (RESULTS.md §2)
that the additive baseline strengthens with the number of contexts averaged, so
a CDI computed over 79 lines is NOT comparable to one computed over 660 --
the larger panel gets a better shared estimate and therefore a smaller
interaction share for purely statistical reasons. Comparing the strata as they
come would manufacture a difference. Every solid-lineage estimate here is
therefore computed on a random subsample of exactly as many lines as the immune
stratum, repeated, and the immune value is compared against that null
distribution. Without this the analysis would be measuring panel size.

Outputs: results/tables/immune_context_prism.csv
         figure bundle results/figures/14_immune/immune_context/
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

from perturbmodel.celldrug import general_sensitivity_by_rep
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PR = ROOT / "data" / "external" / "prism"
FIG = ROOT / "results" / "figures" / "14_immune"
TAB = ROOT / "results" / "tables"
MIN_LINES = 30
N_DRAW = 25
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def decompose_subset(R, K, keep, min_lines):
    """CDI per compound restricted to a given set of cell-line row indices.

    Same replicate-validated estimator used everywhere else in the project:
    leave-one-line-out shared response, interaction as the covariance of
    per-line residuals between independent detection plates at matched dose.

    **General sensitivity is removed first.** A line's response to everything --
    growth rate, seeding density, drug metabolism -- reproduces across disjoint
    compound halves at r = 0.989, and it is shared between replicate detection
    plates just as a genuine interaction is. Left in, it enters the covariance
    and is counted as context-dependence, which inflated every CDI in this
    project by roughly half. The line mean is taken leave-one-compound-out so a
    real interaction cannot be absorbed into it.
    """
    ALPHA, ACPD = general_sensitivity_by_rep(R, K, rows=keep)
    apos = {c: j for j, c in enumerate(ACPD)}

    out = {}
    for cpd, gc in K.groupby("compound", observed=True):
        j = apos[cpd]
        add_n, add_d, cov, npair = 0.0, 0, 0.0, 0
        for dose, gd in gc.groupby("dose", observed=True):
            cols = gd.index.to_numpy()
            if len(cols) < 2:
                continue
            sub = R[np.ix_(keep, cols)]
            with np.errstate(invalid="ignore"):
                lm = np.nanmean(sub, axis=1)
            # each column's correction comes from its OWN plate
            Acol = np.stack([ALPHA[r][:, j] for r in gd.rep.to_numpy()], axis=1)
            valid = np.isfinite(lm) & np.isfinite(Acol).all(axis=1)
            n = int(valid.sum())
            if n < min_lines:
                continue
            loo = (lm[valid].sum() - lm[valid]) / (n - 1)
            add_n += float(np.mean(loo ** 2)) * n; add_d += n
            res = np.full(sub.shape, np.nan, dtype=np.float32)
            res[valid] = sub[valid] - loo[:, None] - Acol[valid]
            for a in range(len(cols)):
                for b in range(a + 1, len(cols)):
                    m = np.isfinite(res[:, a]) & np.isfinite(res[:, b])
                    if m.sum() >= min_lines:
                        cov += float(np.mean(res[m, a] * res[m, b])); npair += 1
        if npair >= 3 and add_d > 0:
            add = add_n / add_d
            inter = max(cov / npair, 0.0)
            den = add + inter
            if den > 0:
                out[cpd] = inter / den
    return pd.Series(out, name="cdi")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-draw", type=int, default=N_DRAW)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

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
    ti["dose_s"] = ti.dose.round(4).astype(str)
    L = lfc.to_numpy(dtype=np.float32)
    colpos = {c: i for i, c in enumerate(lfc.columns)}
    keys, mats = [], []
    for (cpd, dose, rep), g in ti.groupby(["name", "dose_s", "rep"],
                                          observed=True):
        idx = [colpos[c] for c in g.column_name]
        with np.errstate(invalid="ignore"):
            v = np.nanmean(L[:, idx], axis=1) if len(idx) > 1 else L[:, idx[0]]
        keys.append((cpd, dose, rep)); mats.append(v)
    K = pd.DataFrame(keys, columns=["compound", "dose", "rep"])
    R = np.stack(mats, axis=1)

    ci = pd.read_csv(PR / "secondary-screen-cell-line-info.csv")
    lin = dict(zip(ci.depmap_id, ci.primary_tissue.astype(str)))
    sub2 = dict(zip(ci.depmap_id, ci.secondary_tissue.astype(str)))
    tis = np.array([lin.get(d, "unknown") for d in lfc.index])
    imm = np.isin(tis, ["leukemia", "lymphoma"])
    sol = ~imm & ~np.isin(tis, ["unknown", "nan", "fibroblast"])
    print(f"{imm.sum()} immune-lineage lines "
          f"({pd.Series([sub2.get(d,'') for d in lfc.index[imm]]).value_counts().head(4).to_dict()})")
    print(f"{sol.sum()} solid-lineage lines", flush=True)

    n_imm = int(imm.sum())
    cdi_imm = decompose_subset(R, K, np.where(imm)[0], MIN_LINES)
    print(f"\nimmune stratum: {len(cdi_imm)} compounds, median CDI "
          f"{cdi_imm.median():.3f}", flush=True)

    # size-matched solid draws: the null is "a panel of this size", not "solid"
    draws = []
    sol_idx = np.where(sol)[0]
    for d in range(args.n_draw):
        pick = rng.choice(sol_idx, n_imm, replace=False)
        s = decompose_subset(R, K, pick, MIN_LINES)
        draws.append(s.rename(f"d{d}"))
        if d % 5 == 0:
            print(f"  solid draw {d+1}/{args.n_draw}: n={len(s)}, "
                  f"median {s.median():.3f}", flush=True)
    S = pd.concat(draws, axis=1)
    sol_med = S.median(axis=1)
    draw_meds = S.median(axis=0)

    common = cdi_imm.index.intersection(sol_med.index)
    a, b = cdi_imm[common], sol_med[common]
    w = stats.wilcoxon(a, b)
    z = (cdi_imm.median() - draw_meds.mean()) / max(draw_meds.std(), 1e-9)
    print(f"\nsize-matched comparison on {len(common)} shared compounds")
    print(f"  immune median CDI {a.median():.4f}  vs solid (size-matched) "
          f"{b.median():.4f}")
    print(f"  paired Wilcoxon p={w.pvalue:.2e}; immune vs the distribution of "
          f"{args.n_draw} solid draws: z={z:+.2f} "
          f"(solid draws {draw_meds.mean():.4f} +/- {draw_meds.std():.4f})")
    print("  the z uses the spread ACROSS equally sized solid panels, so it "
          "asks whether\n  immune lineage differs by more than panel-to-panel "
          "variation at this n.")

    out = pd.DataFrame({"compound": common, "cdi_immune": a.values,
                        "cdi_solid_matched": b.values})
    moa = ti.drop_duplicates("name").set_index("name").moa
    out["moa"] = out.compound.map(moa)
    out["delta"] = out.cdi_immune - out.cdi_solid_matched
    out.to_csv(TAB / "immune_context_prism.csv", index=False)

    kn = out[out.moa.notna() & (out.moa.astype(str) != "nan")]
    g = (kn.groupby("moa").delta.agg(["median", "size"]).query("size>=5")
         .sort_values("median", ascending=False))
    print(f"\nmechanisms most MORE context-specific in immune lineage "
          f"({len(g)} classes n>=5):")
    print(g.head(8).round(3).to_string())
    print("\n...and most LESS:")
    print(g.tail(6).round(3).to_string())

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    ax[0].hist(draw_meds, bins=12, color="#bdbdbd", label=f"solid panels "
               f"(n={n_imm} lines each)")
    ax[0].axvline(cdi_imm.median(), color=ORANGE, lw=2.5,
                  label="immune lineage")
    ax[0].set_xlabel("median context-dependence index")
    ax[0].set_ylabel("size-matched solid panels")
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title(f"A  Immune vs equally sized solid panels\nz={z:+.2f}",
                    loc="left", fontweight="bold", fontsize=10)

    ax[1].scatter(b, a, s=9, alpha=0.35, color=BLUE, edgecolors="none")
    lim = float(max(a.max(), b.max()))
    ax[1].plot([0, lim], [0, lim], ls="--", color="#555555", lw=1)
    ax[1].set_xlabel("CDI, size-matched solid panel")
    ax[1].set_ylabel("CDI, immune lineage")
    ax[1].set_title(f"B  Per compound (n={len(common)})\n"
                    f"Wilcoxon p={w.pvalue:.1e}", loc="left",
                    fontweight="bold", fontsize=10)

    top = pd.concat([g.head(7), g.tail(7)])
    yy = np.arange(len(top))[::-1]
    ax[2].barh(yy, top["median"],
               color=[ORANGE if v > 0 else AQUA for v in top["median"]],
               height=0.68)
    ax[2].axvline(0, color="#444444", lw=0.9)
    ax[2].set_yticks(yy, [f"{m[:28]} ({int(k)})" for m, k in
                          zip(top.index, top["size"])], fontsize=6.5)
    ax[2].set_xlabel("Δ CDI (immune − size-matched solid)")
    ax[2].set_title("C  Which mechanisms differ by lineage", loc="left",
                    fontweight="bold", fontsize=10)
    fig.suptitle("Is drug response more context-dependent in immune-lineage "
                 "contexts?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "immune_context", FIG,
                    source_data={"per_compound": out,
                                 "by_mechanism": g.reset_index(),
                                 "solid_draw_medians":
                                 draw_meds.rename("median_cdi").reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
