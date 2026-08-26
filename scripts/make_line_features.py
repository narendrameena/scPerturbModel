#!/usr/bin/env python3
"""Phase 3: covariate features for every cell line (for unseen-line context).

Features per line: driver-gene mutation indicators (genes altered in >=3 of the
50 atlas lines) + organ one-hot. Written for ALL 50 atlas lines so held-out-line
models and future dev-subset expansion use the same table.

Output: data/processed/line_features.npz (lines [str], feat [n x d], names [str])
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "processed"

cl = pd.read_parquet(ROOT / "data/metadata/metadata/cell_line_metadata.parquet")
cond = pd.read_csv(ROOT / "results/tables/cells_per_condition.csv")
atlas_lines = sorted(cond.cell_line.unique())
cl = cl[cl.Cell_ID_Cellosaur.isin(atlas_lines)]

gene_counts = cl.groupby("Driver_Gene_Symbol").Cell_ID_Cellosaur.nunique()
genes = sorted(gene_counts[gene_counts >= 3].index)
organs = sorted(cl.drop_duplicates("Cell_ID_Cellosaur").Organ.unique())

feat, names = [], [f"mut:{g}" for g in genes] + [f"organ:{o}" for o in organs]
for ln in atlas_lines:
    rows = cl[cl.Cell_ID_Cellosaur == ln]
    mut = np.isin(genes, rows.Driver_Gene_Symbol.unique()).astype(np.float32)
    org = np.array([o == rows.Organ.iloc[0] for o in organs], np.float32)
    feat.append(np.concatenate([mut, org]))

np.savez_compressed(OUT / "line_features.npz", lines=np.array(atlas_lines),
                    feat=np.stack(feat), names=np.array(names))
print(f"{len(atlas_lines)} lines x {len(names)} features "
      f"({len(genes)} driver genes + {len(organs)} organs) -> line_features.npz")
