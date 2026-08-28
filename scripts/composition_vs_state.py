#!/usr/bin/env python3
"""Is a line x drug interaction a uniform cell-state shift, or a change in the
fraction of cells occupying a PRE-EXISTING state?

This is the question single cells can answer and pseudobulk cannot, and the
field lists it as open ("why do genetically identical cells respond
differently"). Our own results set it up: the interaction is real, reproducible
and multi-dimensional, yet no line-level descriptor predicts it — which is
exactly what one expects if the operative variable is a subpopulation FRACTION
rather than a line average.

For each (line, drug, dose, plate) and each component, treated cells are
compared with plate-matched DMSO cells of the same line:

  uniform shift   every quantile of the distribution moves by the same amount
  compositional   the shift is concentrated in one tail, so quantile shifts
                  have a slope and the variance changes

  shift_q          quantile(treated, q) - quantile(control, q), q = 0.1 .. 0.9
  uniformity       1 - |slope of shift_q across q| * spread / |mean shift|
                   (1 = perfectly uniform, 0 = entirely tail-driven)
  var_ratio        var(treated) / var(control)
  responder_frac   mixture estimate: fraction of treated cells beyond the
                   control distribution's 90th percentile, minus the 10% null
  preexist_frac    fraction of CONTROL cells already beyond the treated
                   population's median -> is the responsive state pre-existing?

Outputs: results/tables/composition_vs_state_<tag>.csv
         figure bundle results/figures/08_single_cell/composition_vs_state_<tag>/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "08_single_cell"
TAB = ROOT / "results" / "tables"
QS = np.arange(0.1, 0.91, 0.1)
MIN_CELLS = 100
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="dev47")
    ap.add_argument("--min-abs-shift", type=float, default=0.02,
                    help="only analyse conditions with a real mean shift")
    args = ap.parse_args()
    FIG.mkdir(parents=True, exist_ok=True)

    src = ROOT / "data" / "processed" / f"cell_component_scores_{args.tag}.parquet"
    df = pd.read_parquet(src)
    comps = [c for c in df.columns if c.startswith("C") and c[1:].isdigit()]
    print(f"{len(df):,} cells, components {comps}", flush=True)

    ctrl_mask = df.drug == "DMSO_TF"
    ctrl_groups = {k: g for k, g in
                   df[ctrl_mask].groupby(["cell_line_id", "plate"],
                                         observed=True)}
    recs = []
    treated = df[~ctrl_mask]
    for (ln, dr, cc, pl), grp in treated.groupby(
            ["cell_line_id", "drug", "conc", "plate"], observed=True):
        c = ctrl_groups.get((ln, pl))
        if c is None or len(grp) < MIN_CELLS or len(c) < MIN_CELLS:
            continue
        for comp in comps:
            t = grp[comp].to_numpy()
            k = c[comp].to_numpy()
            mean_shift = float(t.mean() - k.mean())
            if abs(mean_shift) < args.min_abs_shift:
                continue
            tq = np.quantile(t, QS)
            kq = np.quantile(k, QS)
            shift_q = tq - kq
            # slope of quantile shift: 0 for a pure location shift
            slope = float(np.polyfit(QS, shift_q, 1)[0])
            uniformity = float(1.0 - min(1.0, abs(slope) * 0.8
                                         / (abs(mean_shift) + 1e-9)))
            hi = np.quantile(k, 0.9) if mean_shift > 0 else np.quantile(k, 0.1)
            resp = float((t > hi).mean() - 0.1) if mean_shift > 0 else \
                float((t < hi).mean() - 0.1)
            med_t = float(np.median(t))
            pre = float((k > med_t).mean()) if mean_shift > 0 else \
                float((k < med_t).mean())
            recs.append({
                "cell_line_id": ln, "drug": dr, "conc": cc, "plate": pl,
                "component": comp, "n_treated": len(grp), "n_control": len(c),
                "mean_shift": mean_shift, "quantile_slope": slope,
                "uniformity": uniformity,
                "var_ratio": float(t.var() / max(k.var(), 1e-9)),
                "responder_frac": max(0.0, resp), "preexist_frac": pre,
            })
    res = pd.DataFrame(recs)
    res.to_csv(TAB / f"composition_vs_state_{args.tag}.csv", index=False)
    print(f"\n{len(res):,} (condition, component) comparisons with a shift")
    print(res.groupby("component")[["mean_shift", "uniformity", "var_ratio",
                                    "responder_frac", "preexist_frac"]]
          .median().round(3).to_string())

    frac_comp = float((res.uniformity < 0.5).mean())
    print(f"\ncompositional (uniformity<0.5): {frac_comp:.1%} of comparisons")
    print(f"median responder fraction: {res.responder_frac.median():.3f}")
    print(f"median pre-existing fraction: {res.preexist_frac.median():.3f} "
          f"(0.5 = responsive state fully pre-exists in controls, "
          f"0 = treated cells move outside the control distribution)")

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)

    axes[0].hist(res.uniformity, bins=50, color=BLUE, alpha=0.8)
    axes[0].axvline(0.5, color="#333333", ls="--", lw=1.2)
    axes[0].set_xlabel("uniformity  (1 = all cells shift, 0 = one tail shifts)")
    axes[0].set_ylabel("condition × component")
    axes[0].set_title("A  Uniform shift or subpopulation?", loc="left",
                      fontweight="bold", fontsize=10)
    axes[0].text(0.03, 0.95, f"{frac_comp:.0%} compositional",
                 transform=axes[0].transAxes, fontsize=9, va="top",
                 color="#444444")

    axes[1].scatter(res.responder_frac, res.preexist_frac, s=4, alpha=0.2,
                    color=VIOLET, edgecolors="none")
    axes[1].set_xlabel("responder fraction (treated beyond control tail)")
    axes[1].set_ylabel("pre-existing fraction (control beyond treated median)")
    axes[1].set_title("B  Does the responsive state pre-exist?", loc="left",
                      fontweight="bold", fontsize=10)

    data = [res[res.component == c].uniformity.dropna() for c in comps]
    vp = axes[2].violinplot(data, positions=range(len(comps)), widths=0.7,
                            showmedians=True, showextrema=False)
    for b in vp["bodies"]:
        b.set_facecolor(AQUA); b.set_alpha(0.6); b.set_edgecolor("none")
    vp["cmedians"].set_color("#333333")
    axes[2].axhline(0.5, color="#888888", ls="--", lw=0.9)
    axes[2].set_xticks(range(len(comps)), comps)
    axes[2].set_ylabel("uniformity")
    axes[2].set_title("C  By response program", loc="left", fontweight="bold",
                      fontsize=10)

    fig.suptitle("Is context-specific drug response a uniform cell-state shift "
                 "or a subpopulation effect?", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, f"composition_vs_state_{args.tag}", FIG,
                    source_data=res, script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
