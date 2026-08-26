#!/usr/bin/env python3
"""Expand the dev subset to ALL 47 analyzed cell lines (same 50 drugs).

Motivated by hard-splits result: unseen-line context cannot be learned from
7 training lines. Writes data/interim/dev47/selection.yaml in the same format
consumed by build_dev_subset.py.
"""
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "interim" / "dev47"
OUT.mkdir(parents=True, exist_ok=True)
EXCLUDED = {"NCI-H661", "NCI-H596", "NCI-H2122"}  # paper's low-coverage lines

base = yaml.safe_load(open(ROOT / "data/interim/dev_subset/selection.yaml"))
cond = pd.read_csv(ROOT / "results/tables/cells_per_condition.csv")

lines = (cond[~cond.cell_name.isin(EXCLUDED)]
         .groupby(["cell_line", "cell_name"], observed=True)
         .n_cells.sum().reset_index())
sel = {
    "description": "Expanded dev subset: all 47 analyzed lines x the same 50 "
                   "drugs (+DMSO_TF). Motivation: held-out-line context needs "
                   "more training lines.",
    "cell_lines": [{"cellosaurus": r.cell_line, "name": r.cell_name,
                    "total_cells": int(r.n_cells)}
                   for r in lines.itertuples()],
    "drugs": base["drugs"],
    "controls": ["DMSO_TF"],
}
est = cond[(~cond.cell_name.isin(EXCLUDED))
           & cond.drug.isin(base["drugs"] + ["DMSO_TF"])].n_cells.sum()
sel["estimated_cells"] = int(est)
with open(OUT / "selection.yaml", "w") as fh:
    yaml.safe_dump(sel, fh, sort_keys=False, width=100)
print(f"{len(sel['cell_lines'])} lines x {len(sel['drugs'])} drugs, "
      f"estimated {est:,} cells -> {OUT/'selection.yaml'}")
