#!/usr/bin/env python3
"""Pseudobulk sciPlex3 into (cell_line, perturbation, dose, replicate) profiles.

Restricted to the 24 h timepoint so the comparison with Tahoe-100M (24 h) and
OP3 (24 h) is like-for-like. Groups with fewer than MIN_CELLS cells are dropped,
matching the condition-level QC used throughout.

Output: data/external/scperturb/sciplex3_pseudobulk.h5ad
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data/external/scperturb/sciplex3.h5ad"
OUT = ROOT / "data/external/scperturb/sciplex3_pseudobulk.h5ad"
KEYS = ["cell_line", "perturbation", "dose_value", "replicate"]
MIN_CELLS = 30

A = ad.read_h5ad(SRC)
o = A.obs
keep = (o.time.astype(str) == "24.0") & o.cell_line.astype(str).isin(
    ["MCF7", "A549", "K562"]) & o.replicate.astype(str).isin(["rep1", "rep2"])
A = A[keep.to_numpy()]
o = A.obs
print(f"{A.shape[0]} cells at 24h across {o.cell_line.nunique()} lines, "
      f"{o.perturbation.nunique()} perturbations", flush=True)

X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
grp = o.groupby([o[k].astype(str) for k in KEYS], observed=True).indices
rows, meta = [], []
for k, idx in grp.items():
    if len(idx) < MIN_CELLS:
        continue
    rows.append(np.asarray(X[idx].sum(0)).ravel())
    meta.append(dict(zip(KEYS, k), n_cells=len(idx)))
M = np.stack(rows).astype(np.float32)
obs = pd.DataFrame(meta)
obs.index = obs.index.astype(str)
print(f"{len(obs)} pseudobulk profiles "
      f"({obs.n_cells.median():.0f} median cells)", flush=True)

out = ad.AnnData(X=M, obs=obs, var=A.var.copy())
out.write_h5ad(OUT)
print(f"-> {OUT}")
