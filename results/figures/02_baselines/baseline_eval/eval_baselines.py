#!/usr/bin/env python3
"""Phase 2: pseudobulk baselines under the held-out-condition split.

Setup (log1p CPM space, plate-matched deltas, plate14 excluded):
  profile(L,D,C,plate) - ctrl(L,plate)  ->  cell-weighted mean over plates
  = delta(L,D,C) for every (line, drug, conc) with >=MIN_CELLS on both sides.

Split: random 20% of (L,D,C) triples held out (seed 0); a test triple is kept
only if its (D,C) is still observed in >=1 training line (the additive baseline
needs it; unseen-drug generalization is a separate, harder split).

Baselines:
  no_change      predicted delta = 0 (RMSE only; r undefined for a constant)
  random_delta   delta of a random TRAIN condition of the same line (neg. control)
  additive_shift predicted delta(L,D,C) = mean over train lines of delta(.,D,C)

Metrics per test condition:
  r_hvg   Pearson r of predicted vs true delta over the 2,000 most
          delta-variable genes (variance computed on TRAIN deltas only)
  r_de100 Pearson r over the condition's own top-100 |true delta| genes
  rmse_hvg

Outputs: results/tables/baseline_eval.csv + figure bundle results/figures/baseline_eval/
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from perturbmodel.evaluation import held_out_condition_triples
from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
PB = ROOT / "data" / "processed" / "pseudobulk_dev"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
MIN_CELLS = 50
N_HVG = 2000
TEST_FRAC = 0.2
BLUE, ORANGE = "#2a78d6", "#eb6834"

# ---------------- load & normalize ----------------
counts = np.load(PB / "pseudobulk_counts.npz")["counts"]
cond = pd.read_csv(PB / "conditions.csv")
X = np.log1p(counts / (counts.sum(axis=1, keepdims=True) + 1e-9) * 1e6)

cond = cond.reset_index().rename(columns={"index": "row"})
cond = cond[cond.plate != "plate14"]                       # replicate plate held out
ctrl = cond[(cond.drug == "DMSO_TF") & (cond.n_cells >= MIN_CELLS)]
ctrl_idx = {(r.cell_line_id, r.plate): r.row for r in ctrl.itertuples()}

# ---------------- plate-matched deltas ----------------
rows = []
treated = cond[(cond.drug != "DMSO_TF") & (cond.n_cells >= MIN_CELLS)]
for r in treated.itertuples():
    ci = ctrl_idx.get((r.cell_line_id, r.plate))
    if ci is None:
        continue
    rows.append((r.cell_line_id, r.drug, r.conc, r.n_cells, X[r.row] - X[ci]))
key_df = pd.DataFrame([(a, b, c, d) for a, b, c, d, _ in rows],
                      columns=["cell_line_id", "drug", "conc", "n_cells"])
deltas_all = np.stack([v for *_, v in rows])

# cell-weighted mean across plates per (L,D,C)
key_df["gid"] = key_df.groupby(["cell_line_id", "drug", "conc"], observed=True).ngroup()
n_grp = key_df.gid.nunique()
W = np.zeros((n_grp, deltas_all.shape[1]), dtype=np.float64)
wsum = np.zeros(n_grp)
for i, (gid, w) in enumerate(zip(key_df.gid, key_df.n_cells)):
    W[gid] += w * deltas_all[i]
    wsum[gid] += w
DELTA = (W / wsum[:, None]).astype(np.float32)
G = (key_df.drop_duplicates("gid").sort_values("gid")
     [["cell_line_id", "drug", "conc"]].reset_index(drop=True))
print(f"{len(G)} (line,drug,conc) conditions with plate-matched deltas")

# ---------------- split (shared with single-cell models) ----------------
test_triples = held_out_condition_triples(G, seed=0, test_frac=TEST_FRAC)
test_mask = np.array([(l, d, c) in test_triples
                      for l, d, c in zip(G.cell_line_id, G.drug, G.conc)])
print(f"train {np.sum(~test_mask)}, test {test_mask.sum()} conditions")

hvg = np.argsort(DELTA[~test_mask].var(axis=0))[-N_HVG:]

# additive-shift lookup from TRAIN only (test (L,D,C) never contributes:
# its row is excluded from train by the split itself)
train_G = G[~test_mask].reset_index(drop=True)
train_D = DELTA[~test_mask]
shift = {k: train_D[grp.index.to_numpy()].mean(axis=0)
         for k, grp in train_G.groupby(["drug", "conc"], observed=True)}

# ---------------- evaluate ----------------
def pearson(a, b):
    a, b = a - a.mean(), b - b.mean()
    den = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / den) if den > 0 else np.nan

rng = np.random.default_rng(1)
by_line = {ln: grp.index.to_numpy()
           for ln, grp in train_G.groupby("cell_line_id", observed=True)}

recs = []
for i in np.where(test_mask)[0]:
    true = DELTA[i]
    de100 = np.argsort(np.abs(true))[-100:]
    pool = by_line[G.cell_line_id[i]]
    rand_pred = train_D[rng.choice(pool)]
    for name, pred in (("no_change", np.zeros_like(true)),
                       ("random_delta", rand_pred),
                       ("additive_shift", shift[(G.drug[i], G.conc[i])])):
        recs.append({
            "cell_line_id": G.cell_line_id[i], "drug": G.drug[i], "conc": G.conc[i],
            "baseline": name,
            "r_hvg": pearson(pred[hvg], true[hvg]),
            "r_de100": pearson(pred[de100], true[de100]),
            "rmse_hvg": float(np.sqrt(np.mean((pred[hvg] - true[hvg]) ** 2))),
            "true_norm_hvg": float(np.linalg.norm(true[hvg])),
        })
res = pd.DataFrame(recs)
res.to_csv(TAB / "baseline_eval.csv", index=False)
summary = (res.groupby("baseline")[["r_hvg", "r_de100", "rmse_hvg"]]
           .agg(["mean", "median"]).round(3))
print(summary.to_string())

# ---------------- figure ----------------
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": 0.25, "figure.facecolor": "white"})
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
for ax, metric, title in ((axes[0], "r_hvg", "A  r on 2,000 most responsive genes"),
                          (axes[1], "r_de100", "B  r on condition's top-100 DE genes")):
    data = [res[res.baseline == b][metric].dropna()
            for b in ("random_delta", "additive_shift")]
    vp = ax.violinplot(data, positions=[0, 1], widths=0.7, showmedians=True,
                       showextrema=False)
    for body, c in zip(vp["bodies"], (ORANGE, BLUE)):
        body.set_facecolor(c); body.set_alpha(0.55); body.set_edgecolor("none")
    vp["cmedians"].set_color("#333333")
    ax.axhline(0, color="#888888", lw=0.8, ls="--", zorder=0)
    ax.set_xticks([0, 1], ["random train delta\n(same line)", "additive shift\n(same drug+dose)"])
    ax.set_ylabel("Pearson r (predicted vs true delta)")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=10)
fig.suptitle("Held-out-condition baselines, dev subset pseudobulk",
             fontsize=11, x=0.01, ha="left")
d = save_figure(fig, "baseline_eval", FIG,
                source_data={"baseline_eval": res}, script=__file__)
print(f"figure bundle -> {d}")
