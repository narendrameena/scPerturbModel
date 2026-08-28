#!/usr/bin/env python3
"""Aggregate the dev subset into pseudobulk raw-count profiles.

One profile per (cell_line_id, drug, conc, plate): summed UMI counts over cells
(Squair et al. 2021 style). Outputs:
  data/processed/pseudobulk_dev/pseudobulk_counts.npz   (float32, conds x 62,710 genes)
  data/processed/pseudobulk_dev/conditions.csv          (row index + n_cells)

Usage: python scripts/build_pseudobulk.py [--src data/interim/dev_subset/data]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata" / "metadata"


def token_lut() -> tuple[np.ndarray, pd.DataFrame]:
    gm = pd.read_parquet(META / "gene_metadata.parquet").sort_values("token_id")
    lut = np.full(int(gm.token_id.max()) + 1, -1, dtype=np.int64)
    lut[gm.token_id.to_numpy()] = np.arange(len(gm))
    return lut, gm.reset_index(drop=True)


def dose_lookup() -> dict:
    """sample -> concentration, parsed from sample_metadata (raw mode only)."""
    import ast
    sm = pd.read_parquet(META / "sample_metadata.parquet")
    out = {}
    for s, dnc in zip(sm["sample"], sm.drugname_drugconc):
        try:
            out[s] = float(ast.literal_eval(dnc)[0][1])
        except Exception:
            out[s] = float("nan")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="data/interim/dev_subset/data")
    ap.add_argument("--out", default="data/processed/pseudobulk_dev")
    ap.add_argument("--raw", action="store_true",
                    help="read the untouched atlas shards (all 379 drugs, all "
                         "lines) and join dose from sample_metadata on the fly, "
                         "instead of a pre-extracted subset")
    args = ap.parse_args()
    src, out = ROOT / args.src, ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)

    lut, gm = token_lut()
    n_genes = len(gm)
    acc: dict[tuple, np.ndarray] = {}
    ncells: dict[tuple, int] = {}
    doses = dose_lookup() if args.raw else None

    shards = sorted(src.glob("train-*.parquet"))
    for si, shard in enumerate(shards, 1):
        cols = ["genes", "expressions", "cell_line_id", "drug", "plate"]
        cols += ["sample"] if args.raw else ["conc"]
        t = pq.read_table(shard, columns=cols)
        df = t.to_pandas()
        if args.raw:
            df["conc"] = df["sample"].map(doses).astype("float32")
        for key, grp in df.groupby(["cell_line_id", "drug", "conc", "plate"],
                                   observed=True):
            genes = [np.asarray(g[1:] if e[0] < 0 else g, dtype=np.int64)
                     for g, e in zip(grp.genes, grp.expressions)]
            exprs = [np.asarray(e[1:] if e[0] < 0 else e, dtype=np.float32)
                     for e in grp.expressions]
            G, E = np.concatenate(genes), np.concatenate(exprs)
            cols = lut[G]
            ok = cols >= 0
            vec = acc.setdefault(key, np.zeros(n_genes, dtype=np.float32))
            np.add.at(vec, cols[ok], E[ok])
            ncells[key] = ncells.get(key, 0) + len(grp)
        if si % 50 == 0 or si == len(shards):
            print(f"[{si}/{len(shards)}] {len(acc)} conditions", flush=True)

    keys = sorted(acc)
    mat = np.stack([acc[k] for k in keys]).astype(np.float32)
    cond = pd.DataFrame(keys, columns=["cell_line_id", "drug", "conc", "plate"])
    cond["n_cells"] = [ncells[k] for k in keys]
    np.savez_compressed(out / "pseudobulk_counts.npz", counts=mat)
    cond.to_csv(out / "conditions.csv", index=False)
    gm.to_csv(out / "genes.csv", index=False)
    print(f"pseudobulk: {mat.shape[0]} conditions x {mat.shape[1]} genes "
          f"({mat.nbytes/1e9:.2f} GB) -> {out}")


if __name__ == "__main__":
    main()
