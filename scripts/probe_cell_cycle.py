#!/usr/bin/env python3
"""Phase 1 sanity probe: drug-induced cell-cycle shifts (paper Fig 6C-E),
computed across all 47 analyzed cell lines directly from obs_metadata
(per-cell `phase` is precomputed upstream) — no expression data needed.

Per (cell line, plate, drug, conc) with >=50 treated and >=50 control cells:
  log2 OR(phase) = log2[ (n_phase_t/n_other_t) / (n_phase_c/n_other_c) ]
(0.5 Haldane correction), control = plate-matched DMSO_TF of the same line.

Expected from the paper: palbociclib G1 up; Dinaciclib G2M up; Belinostat/
Panobinostat G2M up; Tucidinostat G2M down; Carbamazepine ~0; microtubule
inhibitors G2M up except Paclitaxel/Tubulin-Inhibitor-6 minimal.

Outputs: results/tables/cell_cycle_log2or.csv, results/figures/cell_cycle_probes.png
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from perturbmodel.utils import save_figure

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata" / "metadata"
FIG = ROOT / "results" / "figures" / "01_paper_replication"
TAB = ROOT / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

EXCLUDED = {"NCI-H661", "NCI-H596", "NCI-H2122"}
MIN_CELLS = 50
PHASES = ["G1", "S", "G2M"]
# validated categorical palette (light mode), slots 1-3
PHASE_COLORS = {"G1": "#2a78d6", "S": "#eb6834", "G2M": "#1baf7a"}
BLUE = "#2a78d6"

# ---------- 1. aggregate phase counts per condition (single pass, columnar) ----------
print("aggregating obs_metadata ...")
pf = pq.ParquetFile(META / "obs_metadata.parquet")
cols = ["cell_name", "plate", "drug", "drugname_drugconc", "phase", "pass_filter"]
parts = []
for batch in pf.iter_batches(columns=cols, batch_size=5_000_000):
    df = batch.to_pandas()
    df = df[(df.pass_filter == "full") & (~df.cell_name.isin(EXCLUDED))]
    parts.append(df.groupby(["cell_name", "plate", "drug", "drugname_drugconc",
                             "phase"], observed=True).size())
counts = (pd.concat(parts).groupby(level=[0, 1, 2, 3, 4], observed=True).sum()
          .rename("n").reset_index())
print(f"aggregated: {len(counts):,} (line,plate,condition,phase) rows")

wide = (counts.pivot_table(index=["cell_name", "plate", "drug", "drugname_drugconc"],
                           columns="phase", values="n", fill_value=0,
                           aggfunc="sum", observed=True)
        .reindex(columns=PHASES, fill_value=0).reset_index())
wide["total"] = wide[PHASES].sum(axis=1)

ctrl = (wide[wide.drug == "DMSO_TF"]
        .groupby(["cell_name", "plate"], observed=True)[PHASES + ["total"]].sum())

treated = wide[(wide.drug != "DMSO_TF") & (wide.total >= MIN_CELLS)].copy()
rows = []
for r in treated.itertuples():
    key = (r.cell_name, r.plate)
    if key not in ctrl.index:
        continue
    c = ctrl.loc[key]
    if c["total"] < MIN_CELLS:
        continue
    for ph in PHASES:
        t_ph, t_ot = getattr(r, ph) + 0.5, (r.total - getattr(r, ph)) + 0.5
        c_ph, c_ot = c[ph] + 0.5, (c["total"] - c[ph]) + 0.5
        rows.append((r.cell_name, r.plate, r.drug, r.drugname_drugconc, ph,
                     np.log2((t_ph / t_ot) / (c_ph / c_ot)), r.total, c["total"]))
lor = pd.DataFrame(rows, columns=["cell_name", "plate", "drug", "drugname_drugconc",
                                  "phase", "log2_or", "n_treated", "n_control"])
dr = pd.read_parquet(META / "drug_metadata.parquet")[["drug", "moa-fine"]]
lor = lor.merge(dr, on="drug", how="left")
lor.to_csv(TAB / "cell_cycle_log2or.csv", index=False)
print(f"log2 OR table: {len(lor):,} rows -> {TAB/'cell_cycle_log2or.csv'}")

# ---------- 2. figure ----------
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": 0.25, "grid.linewidth": 0.5,
                     "figure.facecolor": "white"})

def violins(ax, groups, labels, colors, title, ylabel):
    data = [np.asarray(g) for g in groups]
    pos = np.arange(len(data))
    vp = ax.violinplot(data, positions=pos, widths=0.75, showmedians=True,
                       showextrema=False)
    for body, c in zip(vp["bodies"], colors):
        body.set_facecolor(c); body.set_alpha(0.55); body.set_edgecolor("none")
    vp["cmedians"].set_color("#333333"); vp["cmedians"].set_linewidth(1.2)
    ax.axhline(0, color="#888888", lw=0.8, ls="--", zorder=0)
    ax.set_xticks(pos, labels, rotation=35, ha="right")
    ax.set_title(title, fontsize=10, loc="left", fontweight="bold")
    ax.set_ylabel(ylabel)

fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)

# (A) CDK inhibitors: per drug x phase
cdk_drugs = [d for d in ["palbociclib", "Dinaciclib", "Abemaciclib", "Ribociclib"]
             if d in set(lor.drug)]
groups, labels, colors = [], [], []
for d in cdk_drugs:
    for ph in PHASES:
        v = lor[(lor.drug == d) & (lor.phase == ph)].log2_or
        groups.append(v); labels.append(f"{d[:6]}·{ph}"); colors.append(PHASE_COLORS[ph])
violins(axes[0], groups, labels, colors,
        "A  CDK inhibitors, by phase", "log2 OR vs DMSO")
handles = [plt.Rectangle((0, 0), 1, 1, fc=PHASE_COLORS[p], alpha=0.55) for p in PHASES]
axes[0].legend(handles, PHASES, frameon=False, fontsize=8, loc="upper right")

# (B) HDAC inhibitors: G2M only
hdac = [d for d in ["Belinostat", "Panobinostat", "Tucidinostat", "Carbamazepine"]
        if d in set(lor.drug)]
violins(axes[1], [lor[(lor.drug == d) & (lor.phase == "G2M")].log2_or for d in hdac],
        hdac, [BLUE] * len(hdac), "B  HDAC inhibitors, G2M", "log2 OR (G2M) vs DMSO")

# (C) microtubule inhibitors: G2M only
mts = sorted(lor[lor["moa-fine"] == "Microtubule inhibitor"].drug.unique())
violins(axes[2], [lor[(lor.drug == d) & (lor.phase == "G2M")].log2_or for d in mts],
        mts, [BLUE] * len(mts), "C  Microtubule inhibitors, G2M", "log2 OR (G2M) vs DMSO")

fig.suptitle("Cell-cycle phase shifts vs plate-matched DMSO controls "
             "(47 lines, all doses; replication of Tahoe-100M Fig 6C-E)",
             fontsize=11, x=0.01, ha="left")
# medians read-out + figure bundle (png/svg/pdf + source data + script)
sub = lor[lor.drug.isin(cdk_drugs + hdac + list(mts))]
med = (sub.groupby(["drug", "phase"], observed=True).log2_or.median()
       .unstack()[PHASES].round(2))
print("\nmedian log2 OR per drug:\n", med.to_string())
d = save_figure(fig, "cell_cycle_probes", FIG,
                source_data={"log2or_plotted_drugs": sub,
                             "median_log2or": med.reset_index()},
                script=__file__)
print(f"figure bundle -> {d}")
