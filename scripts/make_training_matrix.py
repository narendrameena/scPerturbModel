#!/usr/bin/env python3
"""Phase 2 prep: convert the dev subset into a training-ready AnnData.

Two vectorized passes over data/interim/dev_subset/data/:
  1. sample cells -> per-gene variance of log1p CPM -> top N_HVG genes
  2. build CSR matrix (cells x HVG) + obs table with the shared
     held-out-condition split column ('train'/'val'/'test'; plate14 rows are
     additionally flagged is_replicate_plate for technical-replication tests)

Output: data/processed/dev_subset_hvg.h5ad (+ hvg_genes.csv)
Run AFTER heavy jobs finish (reads all 3,324 subset shards).
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse

from perturbmodel.evaluation import held_out_condition_triples, split_column

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data" / "interim" / "dev_subset" / "data"
META = ROOT / "data" / "metadata" / "metadata"
OUT = ROOT / "data" / "processed"
OBS_COLS = ["drug", "conc", "conc_unit", "cell_line_id", "plate", "sample",
            "BARCODE_SUB_LIB_ID", "moa-fine"]


def flatten(tbl):
    """Vectorized (genes, expressions) -> flat arrays with CLS stripped."""
    g = tbl["genes"].combine_chunks()
    e = tbl["expressions"].combine_chunks()
    offs = g.offsets.to_numpy()
    gv, ev = g.values.to_numpy(), e.values.to_numpy()
    first = offs[:-1]
    assert (ev[first] < 0).all(), "CLS marker missing in some rows"
    keep = np.ones(len(gv), bool)
    keep[first] = False
    row_len = np.diff(offs) - 1
    return gv[keep], ev[keep], row_len


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-hvg", type=int, default=4000)
    ap.add_argument("--sample-shards", type=int, default=400,
                    help="shards used for HVG estimation")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    gm = pd.read_parquet(META / "gene_metadata.parquet").sort_values("token_id")
    lut = np.full(int(gm.token_id.max()) + 1, -1, dtype=np.int64)
    lut[gm.token_id.to_numpy()] = np.arange(len(gm))
    n_genes = len(gm)
    shards = sorted(SRC.glob("train-*.parquet"))

    # ---------- pass 1: HVG on a shard sample ----------
    rng = np.random.default_rng(0)
    pick = rng.choice(len(shards), min(args.sample_shards, len(shards)),
                      replace=False)
    s1 = np.zeros(n_genes)
    s2 = np.zeros(n_genes)
    n_cells = 0
    for i in pick:
        t = pq.read_table(shards[i], columns=["genes", "expressions"])
        gv, ev, row_len = flatten(t)
        rows = np.repeat(np.arange(len(row_len)), row_len)
        tot = np.bincount(rows, weights=ev, minlength=len(row_len))
        x = np.log1p(ev / tot[rows] * 1e4)          # per-cell CPM(1e4) log1p
        cols = lut[gv]
        np.add.at(s1, cols, x)
        np.add.at(s2, cols, x * x)
        n_cells += len(row_len)
    var = s2 / n_cells - (s1 / n_cells) ** 2         # zeros contribute implicitly
    hvg_cols = np.sort(np.argsort(var)[-args.n_hvg:])
    hvg_genes = gm.iloc[hvg_cols]
    hvg_genes.to_csv(OUT / "hvg_genes.csv", index=False)
    hvg_lut = np.full(n_genes, -1, dtype=np.int64)
    hvg_lut[hvg_cols] = np.arange(len(hvg_cols))
    print(f"HVG: top {args.n_hvg} from {n_cells:,} sampled cells")

    # ---------- pass 2: CSR build ----------
    data_parts, indices_parts, indptr = [], [], [np.array([0], dtype=np.int64)]
    obs_parts = []
    nnz = 0
    for si, shard in enumerate(shards, 1):
        t = pq.read_table(shard)
        gv, ev, row_len = flatten(t)
        cols = hvg_lut[lut[gv]]
        keep = cols >= 0
        rows = np.repeat(np.arange(len(row_len)), row_len)[keep]
        counts = np.bincount(rows, minlength=len(row_len))
        data_parts.append(ev[keep].astype(np.float32))
        indices_parts.append(cols[keep].astype(np.int32))
        indptr.append(nnz + np.cumsum(counts, dtype=np.int64))
        nnz += counts.sum()
        obs_parts.append(t.select(OBS_COLS).to_pandas())
        if si % 300 == 0 or si == len(shards):
            print(f"[{si}/{len(shards)}] nnz={nnz:,}", flush=True)

    X = sparse.csr_matrix(
        (np.concatenate(data_parts), np.concatenate(indices_parts),
         np.concatenate(indptr)),
        shape=(sum(len(o) for o in obs_parts), len(hvg_cols)))
    obs = pd.concat(obs_parts, ignore_index=True)
    obs["is_replicate_plate"] = obs.plate == "plate14"
    test_triples = held_out_condition_triples(obs, seed=0, test_frac=0.2)
    obs["split"] = split_column(obs, test_triples, seed=0, val_frac=0.1)

    adata = ad.AnnData(X=X, obs=obs, var=hvg_genes.set_index("gene_symbol"))
    adata.write_h5ad(OUT / "dev_subset_hvg.h5ad")
    print(f"{adata.shape[0]:,} cells x {adata.shape[1]} HVGs "
          f"({X.data.nbytes/1e9:.1f} GB data) -> {OUT/'dev_subset_hvg.h5ad'}")
    print(obs.split.value_counts().to_string())


if __name__ == "__main__":
    main()
