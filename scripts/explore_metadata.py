#!/usr/bin/env python3
"""Phase 0: summarize Tahoe-100M metadata tables.

Reads the tables downloaded by `download_tahoe100m.py --mode metadata` and writes
summary CSVs to results/tables/. obs_metadata (95.6M rows) is only schema-inspected
plus condition-level aggregation via pyarrow, never fully loaded into pandas.
"""
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "metadata" / "metadata"
OUT = ROOT / "results" / "tables"
OUT.mkdir(parents=True, exist_ok=True)
pd.set_option("display.width", 200, "display.max_columns", 50)


def sec(title):
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


# ---------------- cell lines ----------------
sec("cell_line_metadata")
cl = pd.read_parquet(META / "cell_line_metadata.parquet")
print(f"shape={cl.shape}  columns={list(cl.columns)}")
n_lines = cl["Cell_ID_Cellosaur"].nunique()
print(f"unique cell lines: {n_lines}")
print("\nBy organ (unique lines):")
print(cl.groupby("Organ")["Cell_ID_Cellosaur"].nunique().sort_values(ascending=False))
drivers = (cl.groupby("Driver_Gene_Symbol")["Cell_ID_Cellosaur"].nunique()
           .sort_values(ascending=False))
print("\nTop 15 driver genes (lines carrying an alteration):")
print(drivers.head(15))
cl.to_csv(OUT / "cell_line_metadata.csv", index=False)

# ---------------- drugs ----------------
sec("drug_metadata")
dr = pd.read_parquet(META / "drug_metadata.parquet")
print(f"shape={dr.shape}  columns={list(dr.columns)}")
print(f"unique drugs: {dr['drug'].nunique()}")
print("\nmoa-fine counts:")
print(dr["moa-fine"].value_counts())
print("\nmoa-broad counts:")
print(dr["moa-broad"].value_counts())
print(f"\nhuman-approved: {dict(dr['human-approved'].value_counts())}")
print(f"has SMILES: {(dr['canonical_smiles'].astype(str).str.len() > 0).sum()}")
dr.to_csv(OUT / "drug_metadata.csv", index=False)

# ---------------- samples (doses live here) ----------------
sec("sample_metadata")
sm = pd.read_parquet(META / "sample_metadata.parquet")
print(f"shape={sm.shape}  columns={list(sm.columns)}")
print(sm.head(3))
print(f"\nunique samples: {sm['sample'].nunique()}, unique drugs: {sm['drug'].nunique()}")
print(f"unique (drug, conc) conditions: {sm['drugname_drugconc'].nunique()}")
print("\nsamples per plate:")
print(sm["plate"].value_counts().sort_index())
sm.to_csv(OUT / "sample_metadata.csv", index=False)

# ---------------- genes ----------------
sec("gene_metadata")
gm = pd.read_parquet(META / "gene_metadata.parquet")
print(f"shape={gm.shape}  columns={list(gm.columns)}")
print(gm.head(3))
gm.to_csv(OUT / "gene_metadata.csv", index=False)

# ---------------- obs (per-cell, 95.6M rows) ----------------
sec("obs_metadata (schema + condition-level aggregation)")
pf = pq.ParquetFile(META / "obs_metadata.parquet")
print(f"rows={pf.metadata.num_rows:,}  row_groups={pf.metadata.num_row_groups}")
print("schema:")
for f in pf.schema_arrow:
    print(f"  {f.name}: {f.type}")

# aggregate cells per (cell_line, drug, plate) without loading everything
cols = [c for c in ("cell_line", "cell_name", "cell_line_id", "drug", "plate", "sample")
        if c in pf.schema_arrow.names]
counts = {}
for batch in pf.iter_batches(columns=cols, batch_size=2_000_000):
    df = batch.to_pandas()
    key_cols = [c for c in cols if c != "sample"]
    grp = df.groupby(key_cols, observed=True).size()
    for k, v in grp.items():
        counts[k] = counts.get(k, 0) + v
cond = pd.Series(counts).rename("n_cells").reset_index()
cond.columns = key_cols + ["n_cells"]
cond.to_csv(OUT / "cells_per_condition.csv", index=False)
print(f"\ncondition table -> results/tables/cells_per_condition.csv  ({len(cond):,} rows)")
print(cond["n_cells"].describe())

print("\nDone. Summaries in results/tables/")
