#!/usr/bin/env python3
"""Project every single cell onto the context-response components.

The components were derived from pseudobulk residuals; this scores the
individual cells that make up those pseudobulks, so we can ask whether a
line x drug interaction is a uniform shift of all cells or a change in the
fraction of cells occupying a pre-existing state.

Streaming by design: only 6 numbers per cell are kept, so 20M cells cost
~0.5 GB rather than materialising another expression matrix.

Output: data/processed/cell_component_scores_<tag>.parquet with one row per
cell — component scores plus cell_line_id, drug, conc, plate.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata" / "metadata"
TAB = ROOT / "results" / "tables"
OBS_COLS = ["cell_line_id", "drug", "conc", "plate"]


def build_weights(tag: str, n_comp: int):
    """(n_tokens+1, n_comp) loading matrix indexed by gene token id."""
    load = pd.read_csv(TAB / f"rank1_component_loadings_{tag}.csv")
    stats = pd.read_csv(TAB / f"rank1_viability_{tag}.csv")
    good = stats[stats.split_half_r > 0.2].component.tolist()[:n_comp]
    gm = pd.read_parquet(META / "gene_metadata.parquet").sort_values("token_id")
    # component loadings are ordered like the responsive-gene subset, which was
    # itself taken from the token-sorted gene table; match on gene symbol
    pos = pd.Series(np.arange(len(gm)), index=gm.gene_symbol.to_numpy())
    tok = gm.token_id.to_numpy()
    W = np.zeros((int(tok.max()) + 1, len(good)), dtype=np.float32)
    hit = 0
    for _, r in load.iterrows():
        p = pos.get(r.gene_symbol)
        if p is None or isinstance(p, pd.Series):
            continue
        W[tok[int(p)]] = [r[f"comp{c}"] for c in good]
        hit += 1
    print(f"components {good}; {hit}/{len(load)} genes mapped to tokens")
    return W, good


def score_shard(args):
    shard, W_path, n_comp = args
    W = np.load(W_path)["W"]
    t = pq.read_table(shard, columns=["genes", "expressions"] + OBS_COLS)
    g = t["genes"].combine_chunks()
    e = t["expressions"].combine_chunks()
    offs = g.offsets.to_numpy()
    gv, ev = g.values.to_numpy(), e.values.to_numpy()
    first = offs[:-1]
    keep = np.ones(len(gv), bool)
    keep[first] = False                      # strip the CLS marker token
    gv, ev = gv[keep], ev[keep]
    row_len = np.diff(offs) - 1
    rows = np.repeat(np.arange(len(row_len)), row_len)
    tot = np.bincount(rows, weights=ev, minlength=len(row_len))
    x = np.log1p(ev / np.maximum(tot[rows], 1) * 1e4).astype(np.float32)
    gv = np.clip(gv, 0, W.shape[0] - 1)
    scores = np.zeros((len(row_len), W.shape[1]), dtype=np.float32)
    for k in range(W.shape[1]):
        np.add.at(scores[:, k], rows, x * W[gv, k])
    out = t.select(OBS_COLS).to_pandas()
    for k in range(W.shape[1]):
        out[f"C{k + 1}"] = scores[:, k]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="data/interim/dev47/data")
    ap.add_argument("--tag", default="dev47")
    ap.add_argument("--n-comp", type=int, default=6)
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()

    W, good = build_weights(args.tag, args.n_comp)
    wpath = ROOT / "data" / "processed" / f"_component_weights_{args.tag}.npz"
    np.savez(wpath, W=W)

    shards = sorted((ROOT / args.src).glob("train-*.parquet"))
    print(f"{len(shards)} shards", flush=True)
    parts = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(score_shard, (s, wpath, args.n_comp)) for s in shards]
        for i, f in enumerate(as_completed(futs), 1):
            parts.append(f.result())
            if i % 300 == 0 or i == len(futs):
                print(f"[{i}/{len(futs)}]", flush=True)
    df = pd.concat(parts, ignore_index=True)
    out = ROOT / "data" / "processed" / f"cell_component_scores_{args.tag}.parquet"
    df.to_parquet(out, index=False)
    print(f"{len(df):,} cells scored on {len(good)} components -> {out}")


if __name__ == "__main__":
    main()
