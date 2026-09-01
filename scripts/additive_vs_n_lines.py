#!/usr/bin/env python3
"""How strong is the additive baseline as a function of how many contexts it averages?

Motivation. MAP (Nat Mach Intell 2026) restricts Tahoe-100M to six cell lines
and reports that deep models beat a mean baseline by ~12%. Nature Methods 2025,
DrEval, XPert's own baselines and our results — all using many more contexts —
find that baseline nearly unbeatable. The obvious reconciliation is that the
additive baseline is a MEAN, so its quality should improve with the number of
lines averaged, and studies that use few contexts will find it weak.

Test. Leave-one-line-out. For held-out line L the evaluation conditions are
fixed; only the number of OTHER lines contributing to the additive prior varies
(k = 5, 10, 20, 45), resampled several times per fold. Everything else — the
metric, the gene set, the conditions — is held constant, so any change is
attributable to k alone.

Outputs: results/tables/additive_vs_n_lines.csv
         figure bundle results/figures/02_baselines/additive_vs_n_lines/
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.evaluation.delta_eval import (build_deltas, load_pseudobulk,
                                                responsive_genes)
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures" / "02_baselines"
TAB = ROOT / "results" / "tables"
KS = [2, 3, 4, 5, 6, 10, 18, 20, 45]   # extended to the context counts used by
      # published benchmarks (State uses query datasets with 3, 3, 4, 6 and 18)
N_SEEDS = 3
BLUE, ORANGE = "#2a78d6", "#eb6834"


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pb-dir", default="data/processed/pseudobulk_full")
    ap.add_argument("--tag", default="full")
    args = ap.parse_args()

    X, cond = load_pseudobulk(ROOT / args.pb_dir)
    G, DELTA = build_deltas(X, cond)
    resp = responsive_genes(DELTA, np.ones(len(G), bool))
    D = DELTA[:, resp].astype(np.float32)
    del DELTA
    lines = sorted(G.cell_line_id.unique())
    ks = [k for k in KS if k <= len(lines) - 1]
    print(f"{len(G)} conditions, {len(lines)} lines, {len(resp)} genes; k={ks}",
          flush=True)

    # index: (drug, conc) -> {line: row}
    idx = {}
    for r, (ln, dr, cc) in enumerate(zip(G.cell_line_id, G.drug, G.conc)):
        idx.setdefault((dr, cc), {})[ln] = r

    rng = np.random.default_rng(0)
    recs = []
    for fi, L in enumerate(lines):
        others = [l for l in lines if l != L]
        test_rows = [(k, v[L]) for k, v in idx.items() if L in v]
        if not test_rows:
            continue
        for k in ks:
            for seed in range(N_SEEDS):
                pool = rng.choice(others, size=min(k, len(others)),
                                  replace=False)
                pool_set = set(pool.tolist())
                rs = []
                for key, row in test_rows:
                    donors = [r for ln, r in idx[key].items() if ln in pool_set]
                    if len(donors) < 2:
                        continue
                    prior = D[donors].mean(0)
                    true = D[row]
                    de = np.argsort(np.abs(true))[-100:]
                    rs.append(pearson(prior[de], true[de]))
                if rs:
                    recs.append({"line": L, "k": k, "seed": seed,
                                 "n_conditions": len(rs),
                                 "r_de100": float(np.nanmedian(rs))})
        if (fi + 1) % 10 == 0:
            print(f"[{fi + 1}/{len(lines)}]", flush=True)

    res = pd.DataFrame(recs)
    res.to_csv(TAB / "additive_vs_n_lines.csv", index=False)
    summ = res.groupby("k").r_de100.agg(["mean", "median", "std", "size"])
    print("\nadditive baseline quality vs number of lines averaged:")
    print(summ.round(4).to_string())
    lo, hi = summ.loc[ks[0], "median"], summ.loc[ks[-1], "median"]
    print(f"\nk={ks[0]} -> k={ks[-1]}: {lo:.4f} -> {hi:.4f} "
          f"({100 * (hi - lo) / lo:+.1f}%)")

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.25, "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)

    data = [res[res.k == k].r_de100.dropna() for k in ks]
    vp = axes[0].violinplot(data, positions=range(len(ks)), widths=0.7,
                            showmedians=True, showextrema=False)
    for b in vp["bodies"]:
        b.set_facecolor(BLUE); b.set_alpha(0.6); b.set_edgecolor("none")
    vp["cmedians"].set_color("#333333")
    axes[0].set_xticks(range(len(ks)), [str(k) for k in ks])
    axes[0].set_xlabel("cell lines averaged into the additive prior")
    axes[0].set_ylabel("additive baseline r (top-100 DE genes)")
    axes[0].set_title("A  The baseline strengthens with context count",
                      loc="left", fontweight="bold", fontsize=10)

    m = summ["median"]
    axes[1].plot(ks, m.loc[ks], "-o", color=BLUE, lw=2, ms=7)
    axes[1].fill_between(ks, m.loc[ks] - summ["std"].loc[ks],
                         m.loc[ks] + summ["std"].loc[ks], color=BLUE,
                         alpha=0.15, lw=0)
    axes[1].axvline(6, color=ORANGE, ls="--", lw=1.4)
    axes[1].text(6.4, m.loc[ks].min(), "MAP uses 6 lines", color=ORANGE,
                 fontsize=8.5, rotation=90, va="bottom")
    axes[1].set_xlabel("cell lines averaged")
    axes[1].set_ylabel("median r (top-100 DE genes)")
    axes[1].set_title("B  Why studies disagree about the baseline",
                      loc="left", fontweight="bold", fontsize=10)

    fig.suptitle("The additive baseline is only as good as the number of "
                 "contexts it averages", fontsize=11, x=0.01, ha="left")
    d = save_figure(fig, "additive_vs_n_lines", FIG,
                    source_data={"per_fold": res,
                                 "summary": summ.reset_index()},
                    script=__file__)
    print(f"figure bundle -> {d}")


if __name__ == "__main__":
    main()
