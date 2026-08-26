#!/usr/bin/env python3
"""Phase 1: choose the dev subset (cell lines + drugs) for fast iteration.

Data-driven selection:
  - 8 cell lines: 2 KRAS-G12C, 2 KRAS non-G12C, 2 BRAF-V600E (RAS-wt), 2 RAS/RAF-wt,
    preferring high per-condition coverage and organ diversity.
  - ~50 drugs: the paper's biological-probe drugs (RAS/RAF, CDK, HDAC, microtubule,
    proteasome, translation) + >=1 drug per fine-MOA class + a few 'unclear' MOA
    drugs, preferring coverage. DMSO_TF is always included.

Writes data/interim/dev_subset/selection.yaml with the lists + estimated cell count.
Requires: scripts/explore_metadata.py has produced results/tables/cells_per_condition.csv
"""
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata" / "metadata"
OUT = ROOT / "data" / "interim" / "dev_subset"
OUT.mkdir(parents=True, exist_ok=True)

N_LINES_PER_GROUP = {"KRAS_G12C": 2, "KRAS_other": 2, "BRAF_V600E": 2, "RAS_RAF_WT": 2}
N_DRUGS_TOTAL = 50
N_UNCLEAR_MOA = 5
RAS_RAF_GENES = {"KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "ARAF"}

# Probe drugs anchoring the paper's validated context-dependent effects (Fig 6).
PROBE_DRUGS = [
    "Dabrafenib", "RMC-6236", "Adagrasib",                      # RAS/RAF context
    "palbociclib", "Dinaciclib", "Abemaciclib", "Ribociclib",   # CDK / cell cycle
    "Belinostat", "Panobinostat", "Tucidinostat", "Carbamazepine",  # HDAC spectrum
    "Paclitaxel", "Vinorelbine",                                # microtubule
    "Bortezomib",                                               # proteasome (large effect)
    "Homoharringtonine", "Harringtonine",                       # E-distance outliers
]

cond = pd.read_csv(ROOT / "results" / "tables" / "cells_per_condition.csv")
cl = pd.read_parquet(META / "cell_line_metadata.parquet")
dr = pd.read_parquet(META / "drug_metadata.parquet")

# ---------- cell line groups ----------
lines = (cond.groupby(["cell_line", "cell_name"], observed=True)
         .agg(total_cells=("n_cells", "sum"), n_conds=("n_cells", "size"),
              med_cells=("n_cells", "median")).reset_index())
cl_in = cl[cl.Cell_ID_Cellosaur.isin(lines.cell_line)]

def group_of(cvcl: str) -> str:
    rows = cl_in[cl_in.Cell_ID_Cellosaur == cvcl]
    kras = rows[rows.Driver_Gene_Symbol == "KRAS"]
    braf_v600e = ((rows.Driver_Gene_Symbol == "BRAF") &
                  (rows.Driver_ProtEffect_or_CdnaEffect == "p.V600E")).any()
    if (kras.Driver_ProtEffect_or_CdnaEffect == "p.G12C").any():
        return "KRAS_G12C"
    if len(kras):
        return "KRAS_other"
    if braf_v600e:
        return "BRAF_V600E"
    if not rows.Driver_Gene_Symbol.isin(RAS_RAF_GENES).any():
        return "RAS_RAF_WT"
    return "other_RAS_RAF"  # non-KRAS RAS mutants etc. — skip for clean contrasts

lines["group"] = lines.cell_line.map(group_of)
lines["organ"] = lines.cell_line.map(
    cl_in.drop_duplicates("Cell_ID_Cellosaur").set_index("Cell_ID_Cellosaur").Organ)

picked_lines = []
for grp, n in N_LINES_PER_GROUP.items():
    pool = (lines[lines.group == grp]
            .sort_values(["med_cells", "total_cells"], ascending=False))
    chosen, organs = [], set()
    for _, r in pool.iterrows():          # greedy: coverage first, new organ preferred
        if len(chosen) >= n:
            break
        if r.organ in organs and len(pool) - len(chosen) > 1:
            continue
        chosen.append(r); organs.add(r.organ)
    while len(chosen) < n:                # fallback if organ constraint starved us
        rest = pool[~pool.cell_line.isin([c.cell_line for c in chosen])]
        chosen.append(rest.iloc[0])
    picked_lines += chosen
sel_lines = pd.DataFrame(picked_lines)

# ---------- drugs ----------
drug_cells = (cond[cond.cell_line.isin(sel_lines.cell_line)]
              .groupby("drug", observed=True).n_cells.sum())
dr = dr.assign(cells_in_sel=dr.drug.map(drug_cells).fillna(0))

sel_drugs = []
for p in PROBE_DRUGS:
    hit = dr[dr.drug.str.strip().str.lower() == p.strip().lower()]
    if len(hit):
        sel_drugs.append(hit.iloc[0].drug)
    else:
        print(f"WARNING: probe drug not found: {p}")

for moa, grp in dr[dr["moa-fine"] != "unclear"].groupby("moa-fine"):
    if not grp.drug.isin(sel_drugs).any():
        sel_drugs.append(grp.sort_values("cells_in_sel", ascending=False).iloc[0].drug)

n_left = N_DRUGS_TOTAL - len(sel_drugs)
extra = (dr[(dr["moa-fine"] == "unclear") & ~dr.drug.isin(sel_drugs)]
         .sort_values("cells_in_sel", ascending=False)
         .head(max(N_UNCLEAR_MOA, n_left)))
sel_drugs += extra.drug.tolist()[:max(N_UNCLEAR_MOA, n_left)]

# ---------- estimate + write ----------
mask = (cond.cell_line.isin(sel_lines.cell_line)
        & cond.drug.isin(sel_drugs + ["DMSO_TF"]))
est_cells = int(cond[mask].n_cells.sum())

selection = {
    "description": "Dev subset for fast iteration (Phase 1). All doses of each drug "
                   "+ plate-matched DMSO_TF controls for the selected lines.",
    "estimated_cells": est_cells,
    "cell_lines": [
        {"cellosaurus": r.cell_line, "name": r.cell_name, "group": r.group,
         "organ": r.organ, "total_cells": int(r.total_cells)}
        for r in sel_lines.itertuples()
    ],
    "drugs": sorted(sel_drugs),
    "controls": ["DMSO_TF"],
}
with open(OUT / "selection.yaml", "w") as fh:
    yaml.safe_dump(selection, fh, sort_keys=False, width=100)

print(f"Selected {len(sel_lines)} lines / {len(sel_drugs)} drugs "
      f"(+DMSO_TF), estimated {est_cells:,} cells")
print(sel_lines[["cell_name", "group", "organ", "total_cells", "med_cells"]]
      .to_string(index=False))
moa_cov = dr[dr.drug.isin(sel_drugs)]["moa-fine"].value_counts()
print(f"\nMOA classes covered: {moa_cov.index.nunique()} "
      f"(of {dr['moa-fine'].nunique()})")
print(f"-> {OUT/'selection.yaml'}")
