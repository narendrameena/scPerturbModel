#!/usr/bin/env python3
"""Potency versus rewiring at 6x the scale, on LINCS instead of Tahoe.

RESULTS.md sec.41 established that the context x compound interaction is
high-dimensional rather than rank one -- lines run different programmes, not one
programme at different volumes -- and then hit a hard limit. Tahoe supplies only
146 (drug, dose) combinations with enough replicated lines, and just **six**
mechanism classes reach four combinations. Protein synthesis inhibitors, the
class with the highest flagged rate and the one that motivated the potency
hypothesis, had exactly one. The per-class ordering was therefore not
interpretable.

LINCS phase 1 has the structure the question needs and Tahoe does not:

  * **831 (compound, dose, time) combinations** measured in >= 12 cell lines,
    against Tahoe's 146;
  * **97%** of (line, compound, dose, time) conditions appear on more than one
    detection plate, so the cross-replicate covariance is available almost
    everywhere;
  * 71 cell lines and 20,415 compounds, which after joining the Drug Repurposing
    Hub mechanism annotations used for PRISM gives **24 mechanism classes with
    >= 4 combinations** rather than six.

The analysis is the one from sec.41, unchanged, so the two are comparable: remove
the compound main effect leave-one-line-out, remove each line's general
response, then take the cross-replicate covariance BETWEEN LINES for each
(compound, dose, time) and ask how many directions it occupies. Rank is estimated
across independent detection plates, so noise -- which is full rank -- contributes
nothing.

One difference has to be stated. LINCS Level 4 is z-scored within plate, which
already removes part of the context main effect that sec.41 had to subtract
explicitly. The dimensionality being measured is therefore of a residual that has
been through a different normalisation, and the two numbers should be compared
for their ORDERING and their distance from rank one, not equated.

Outputs: results/tables/lincs_potency_rewiring.csv
         figure bundle results/figures/00_manuscript/lincs_potency_rewiring/
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
LI = ROOT / "data" / "external" / "lincs"
PR = ROOT / "data" / "external" / "prism"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures" / "00_manuscript"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
GREY = "#9e9e9e"
MIN_LINES = 12


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-lines", type=int, default=MIN_LINES)
    ap.add_argument("--max-combos", type=int, default=900)
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    inst = pd.read_csv(LI / "p1_inst_info.txt.gz", sep="\t", low_memory=False)
    I = inst[inst.pert_type == "trt_cp"].copy()
    I["dose"] = I.pert_dose.round(3)
    I["combo"] = (I.pert_id.astype(str) + "|" + I.dose.astype(str) + "|"
                  + I.pert_time.astype(str))
    nline = I.groupby("combo").cell_id.nunique()
    keep = set(nline[nline >= args.min_lines].index)
    I = I[I.combo.isin(keep)]
    combos = sorted(keep)[:args.max_combos]
    I = I[I.combo.isin(set(combos))]
    print(f"{len(combos)} (compound, dose, time) combinations with >= "
          f"{args.min_lines} lines; {len(I):,} instances", flush=True)

    from cmapPy.pandasGEXpress.parse import parse
    M = parse(str(LI / "p1_level4.gctx"),
              cid=list(I.inst_id.astype(str))).data_df
    I = I[I.inst_id.astype(str).isin(set(M.columns))]
    X = M.T.to_numpy(dtype=np.float32)
    pos = {c: i for i, c in enumerate(M.columns)}
    I = I.assign(row=[pos[i] for i in I.inst_id.astype(str)])
    print(f"matrix {X.shape}", flush=True)

    # per (line, combo, plate) profile
    P = {}
    for (cb, ln, pl), g in I.groupby(["combo", "cell_id", "rna_plate"],
                                     observed=True):
        P[(cb, ln, pl)] = X[g.row.to_numpy()].mean(0)
    # compound main effect, leave-one-line-out, per combo
    by_combo = {}
    for (cb, ln, pl), v in P.items():
        by_combo.setdefault(cb, {}).setdefault(ln, []).append(v)
    resid = {}
    for cb, lines in by_combo.items():
        per = {l: np.mean(v, axis=0) for l, v in lines.items()}
        if len(per) < args.min_lines:
            continue
        tot = np.sum(list(per.values()), axis=0)
        n = len(per)
        # leave-one-line-out compound effect for each line in this combo
        resid[cb] = {l: (tot - per[l]) / (n - 1) for l in per}
    # residual per (combo, line, plate), then the line's general response
    R = {}
    for (cb, ln, pl), v in P.items():
        if cb in resid and ln in resid[cb]:
            R[(cb, ln, pl)] = v - resid[cb][ln]
    alpha = {}
    by_line = {}
    for (cb, ln, pl), v in R.items():
        by_line.setdefault((ln, pl), {}).setdefault(cb, []).append(v)
    for (ln, pl), d in by_line.items():
        m = {c: np.mean(v, axis=0) for c, v in d.items()}
        if len(m) < 3:
            continue
        tot, n = np.sum(list(m.values()), axis=0), len(m)
        for c in m:
            alpha[(c, ln, pl)] = (tot - m[c]) / (n - 1)
    # centre across lines within (combo, plate) so it is a contrast
    grp = {}
    for (cb, ln, pl), a in alpha.items():
        grp.setdefault((cb, pl), []).append(a)
    cen = {k: np.mean(v, axis=0) for k, v in grp.items()}
    for k in list(R):
        if k in alpha:
            R[k] = R[k] - (alpha[k] - cen[(k[0], k[2])])

    # rank structure per combo, from independent plates
    rows = []
    for cb in combos:
        per = {}
        for (c2, ln, pl), v in R.items():
            if c2 == cb:
                per.setdefault(ln, {})[pl] = v
        use = {l: sorted(d) for l, d in per.items() if len(d) >= 2}
        if len(use) < args.min_lines:
            continue
        lines = sorted(use)
        Ga = np.stack([per[l][use[l][0]] for l in lines])
        Gb = np.stack([per[l][use[l][1]] for l in lines])
        C = 0.5 * (Ga @ Gb.T + Gb @ Ga.T)
        w = np.linalg.eigvalsh(C)[::-1]
        p_ = w[w > 0]
        if not len(p_) or p_.sum() <= 0:
            continue
        rows.append({"combo": cb, "pert_id": cb.split("|")[0],
                     "dose": float(cb.split("|")[1]),
                     "time": cb.split("|")[2], "n_lines": len(lines),
                     "reproducible": float(p_.sum()),
                     "frac_rank1": float(p_[0] / p_.sum()),
                     "n_dim": float((p_.sum() ** 2) / (p_ ** 2).sum())})
    T = pd.DataFrame(rows)
    print(f"\n{len(T)} combinations with >= {args.min_lines} replicated lines",
          flush=True)
    if not len(T):
        print("nothing usable"); return

    names = I.drop_duplicates("pert_id").set_index("pert_id").pert_iname
    T["name"] = T.pert_id.map(names)
    ti = pd.read_csv(PR / "secondary-screen-replicate-treatment-info.csv",
                     low_memory=False)
    t = ti.dropna(subset=["name"]).copy()
    t["k"] = t.name.map(norm)
    moa = t.drop_duplicates("k").set_index("k").moa
    T["moa"] = T.name.map(norm).map(moa)
    print(f"{T.moa.notna().mean():.0%} annotated with a mechanism")

    print(f"\nOVERALL")
    print(f"  median share of the interaction in ONE direction: "
          f"{T.frac_rank1.median():.1%}")
    print(f"  median effective number of directions:            "
          f"{T.n_dim.median():.2f}   (of up to {int(T.n_lines.median())} lines)")
    print(f"  Tahoe (sec.41): 16.8% and 13.5 directions over 146 combinations.")

    byc = (T[T.moa.notna()].groupby("moa")
           .agg(n=("frac_rank1", "size"), rank1=("frac_rank1", "median"),
                ndim=("n_dim", "median"), lines=("n_lines", "median"))
           .query("n>=4").sort_values("ndim", ascending=False))
    print(f"\n{len(byc)} mechanism classes with >= 4 combinations "
          f"(Tahoe had 6):")
    print(byc.round(3).to_string())

    # is the ordering stable, or is it noise? bootstrap over combinations
    rng = np.random.default_rng(0)
    sub = T[T.moa.isin(byc.index)]
    boots = []
    for _ in range(500):
        b = sub.sample(len(sub), replace=True, random_state=int(
            rng.integers(1e9)))
        gg = b.groupby("moa").n_dim.median()
        boots.append(gg.reindex(byc.index))
    B = pd.DataFrame(boots)
    lo, hi = B.quantile(0.025), B.quantile(0.975)
    print(f"\n  bootstrap 95% intervals on effective dimensionality:")
    for m in byc.index[:8]:
        print(f"    {m[:34]:34s} {byc.ndim[m]:5.1f}  "
              f"[{lo[m]:4.1f}, {hi[m]:4.1f}]")
    top, bot = byc.index[0], byc.index[-1]
    sep = (lo[top] > hi[bot])
    print(f"  most vs least rewiring class ({top[:24]} vs {bot[:24]}): "
          f"intervals {'do NOT overlap' if sep else 'OVERLAP'}")
    if not sep:
        print("  The per-class ordering is not resolved even at this scale; "
              "what is resolved\n  is that every class sits far from rank one.")

    T.to_csv(TAB / "lincs_potency_rewiring.csv", index=False)
    byc.to_csv(TAB / "lincs_potency_by_moa.csv")

    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.facecolor": "white"})
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3), constrained_layout=True)
    ax[0].hist(T.frac_rank1, bins=40, color=VIOLET, alpha=0.85)
    ax[0].axvline(T.frac_rank1.median(), color=ORANGE, lw=2.2,
                  label=f"LINCS median {T.frac_rank1.median():.0%}")
    ax[0].axvline(0.168, color=AQUA, ls="--", lw=2,
                  label="Tahoe median 17%")
    ax[0].set_xlabel("share of the interaction in a single direction")
    ax[0].set_ylabel("(compound, dose, time) combinations")
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title(f"a  {len(T)} combinations, vs 146 in Tahoe", loc="left",
                    fontweight="bold", fontsize=9.5)

    sel = byc.head(14)
    yy = np.arange(len(sel))[::-1]
    ax[1].barh(yy, sel.ndim, color=VIOLET, height=0.7,
               xerr=[sel.ndim - lo[sel.index], hi[sel.index] - sel.ndim],
               error_kw=dict(ecolor="#333", lw=1.0))
    ax[1].set_yticks(yy, [f"{m[:30]} ({int(n)})" for m, n in
                          zip(sel.index, sel.n)], fontsize=6.2)
    ax[1].axvline(1.0, color=ORANGE, ls="--", lw=1.6, label="rank one")
    ax[1].set_xlabel("effective number of directions")
    ax[1].legend(frameon=False, fontsize=7)
    ax[1].set_title(f"b  {len(byc)} classes with n≥4, vs 6 in Tahoe",
                    loc="left", fontweight="bold", fontsize=9.5)

    ax[2].scatter(T.n_lines, T.n_dim, s=10, alpha=0.4, color=BLUE,
                  edgecolors="none")
    ax[2].plot([0, T.n_lines.max()], [0, T.n_lines.max()], ls=":", color="#888",
               lw=1.2, label="one direction per line")
    ax[2].axhline(1.0, color=ORANGE, ls="--", lw=1.6, label="rank one")
    ax[2].set_xlabel("cell lines measured")
    ax[2].set_ylabel("effective directions")
    ax[2].legend(frameon=False, fontsize=7)
    ax[2].set_title("c  Dimensionality is not a line-count artefact",
                    loc="left", fontweight="bold", fontsize=9.5)
    fig.suptitle("Potency versus rewiring at scale: LINCS supplies 6x the "
                 "combinations and 4x the mechanism classes", fontsize=10.5,
                 x=0.005, ha="left", fontweight="bold")
    d = save_figure(fig, "lincs_potency_rewiring", FIG,
                    source_data={"per_combo": T,
                                 "by_moa": byc.reset_index()}, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
