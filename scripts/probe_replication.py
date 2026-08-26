#!/usr/bin/env python3
"""Phase 1 sanity probe: biological-replicate agreement (paper Fig 3D) on the dev
subset. Plate 14 replicates plate 6; matched (line, drug, conc) pseudobulk
profiles should correlate ~0.97, unmatched pairs ~0.91.

Requires: scripts/build_pseudobulk.py output in data/processed/pseudobulk_dev/.
Outputs: results/tables/replication_correlations.csv,
         results/figures/replication_plate6_vs_14.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PB = ROOT / "data" / "processed" / "pseudobulk_dev"
FIG = ROOT / "results" / "figures" / "01_paper_replication"
TAB = ROOT / "results" / "tables"
MIN_CELLS = 50
RNG = np.random.default_rng(0)
BLUE, ORANGE = "#2a78d6", "#eb6834"  # validated categorical slots 1-2

counts = np.load(PB / "pseudobulk_counts.npz")["counts"]
cond = pd.read_csv(PB / "conditions.csv")
lognorm = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e6)

six = cond[(cond.plate == "plate6") & (cond.n_cells >= MIN_CELLS)]
fourteen = cond[(cond.plate == "plate14") & (cond.n_cells >= MIN_CELLS)]
key = ["cell_line_id", "drug", "conc"]
merged = six.reset_index().merge(fourteen.reset_index(), on=key,
                                 suffixes=("_6", "_14"))

def pearson(i, j):
    a, b = lognorm[i], lognorm[j]
    am, bm = a - a.mean(), b - b.mean()
    return float((am @ bm) / (np.linalg.norm(am) * np.linalg.norm(bm) + 1e-12))

matched = np.array([pearson(r.index_6, r.index_14) for r in merged.itertuples()])

n_unmatched = min(4000, 4 * len(matched))
ii = RNG.integers(0, len(six), n_unmatched)
jj = RNG.integers(0, len(fourteen), n_unmatched)
unmatched = []
six_idx, ftn_idx = six.reset_index(), fourteen.reset_index()
for i, j in zip(ii, jj):
    a, b = six_idx.iloc[i], ftn_idx.iloc[j]
    if (a.cell_line_id, a.drug, a.conc) == (b.cell_line_id, b.drug, b.conc):
        continue
    unmatched.append(pearson(a["index"], b["index"]))
unmatched = np.array(unmatched)

q = lambda x: np.percentile(x, [25, 50, 75]).round(3)
print(f"matched   (n={len(matched)}):  q25/med/q75 = {q(matched)}")
print(f"unmatched (n={len(unmatched)}): q25/med/q75 = {q(unmatched)}")

corr_df = pd.DataFrame({
    "kind": ["matched"] * len(matched) + ["unmatched"] * len(unmatched),
    "pearson": np.concatenate([matched, unmatched]),
})
corr_df.to_csv(TAB / "replication_correlations.csv", index=False)

plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": 0.25, "figure.facecolor": "white"})
fig, ax = plt.subplots(figsize=(5.5, 4), constrained_layout=True)
bins = np.linspace(min(unmatched.min(), matched.min()), 1.0, 60)
ax.hist(unmatched, bins=bins, density=True, alpha=0.55, color=ORANGE,
        label=f"unmatched (med {np.median(unmatched):.3f})")
ax.hist(matched, bins=bins, density=True, alpha=0.55, color=BLUE,
        label=f"matched (med {np.median(matched):.3f})")
ax.set_xlabel("Pearson r of pseudobulk log1p CPM (plate 6 vs plate 14)")
ax.set_ylabel("density")
ax.set_title("Replicate-plate agreement, dev subset (cf. paper Fig 3D)",
             loc="left", fontweight="bold", fontsize=10)
ax.legend(frameon=False, fontsize=9)
d = save_figure(fig, "replication_plate6_vs_14", FIG,
                source_data=corr_df, script=__file__)
print(f"figure bundle -> {d}")
