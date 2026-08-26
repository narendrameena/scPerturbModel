# Verified data facts — Tahoe-100M (checked 2026-08-26 against the HF release)

Source: `tahoebio/Tahoe-100M` on Hugging Face. License: **CC0-1.0** (public domain).
Local copy: `data/raw/tahoe-100m/` (429 GB total: 338 GB expression, 89 GB pseudobulk
DE, 2.3 GB core metadata). Core metadata also in `data/metadata/metadata/`.

## Expression data (`data/train-*.parquet`, 3,388 shards)

- **95,624,334 records = exactly the paper's full-filter cells.** Per-cell QC is
  already applied upstream; only condition-level filtering (≥50 cells) is on us.
- Fields per record: `genes` (int64 token IDs), `expressions` (float32 raw counts),
  `drug`, `sample`, `BARCODE_SUB_LIB_ID` (unique cell key), `cell_line_id`
  (Cellosaurus), `moa-fine`, `canonical_smiles`, `pubchem_cid`, `plate`.
- **Quirk:** first entry of `genes`/`expressions` is a CLS/marker token
  (`expressions[0] < 0`, typically -2.0; `genes[0] = 1`) — strip before use
  (`perturbmodel.data.dataset.strip_marker`).
- **Shards are organized by plate** (e.g. shard 0 = plate4 only, ~28k cells/shard).
  Training loaders must shuffle across shards, not just within.
- **Dose is NOT in the expression records** — join `sample` →
  `sample_metadata.drugname_drugconc`, e.g. `[('Infigratinib', 0.05, 'uM')]`.

## Metadata tables (`metadata/`)

| table | rows | notes |
|---|---|---|
| `gene_metadata` | 62,710 | token_id ↔ ensembl_id ↔ gene_symbol (Ensembl 109/GRCh38) |
| `sample_metadata` | 1,344 (= 14 plates × 96 wells) | per-well QC means + **drug concentration**; 1,138 unique (drug, conc) |
| `drug_metadata` | 379 drugs | targets/MOA are **GPT-4o-derived → noisy**; 199/379 have moa-fine "unclear"; all have SMILES + PubChem CID |
| `cell_line_metadata` | 1,000 rows, **102 lines** (long format, one row per driver mutation) | covers more lines than the 50 in the atlas — intersect with obs |
| `obs_metadata` | 100,648,790 cells | includes minimal-filter cells; `pass_filter`: 95,624,334 "full" + 5,024,456 "minimal". Has per-cell `phase`, `S_score`, `G2M_score`, QC counts |
| `pseudobulk_differential_expression` | 1,026 shards, 89 GB | precomputed DE — useful for evaluation without recomputing |

- Controls: `drug == "DMSO_TF"`, present on every plate → plate-matched controls.
- 50 cell lines in obs; 380 drug values (379 drugs + DMSO_TF).
- Paper's excluded low-coverage lines (NCI-H661, NCI-H596, NCI-H2122) ARE present in
  the release — drop them for analyses mirroring the paper.

## Compute environment

- Cluster is **CPU-only** (no GPU partitions in SLURM: partitions 2204/2004/1804/NXFL,
  128-core / 385 GB nodes). Implication: dev-subset-scale models (scVI, conditional
  decoders, latent shift) are fine; transformer-scale pretraining would need outside
  GPUs. Env is a **venv** at `.venv/` (no conda on cluster), shared via BeeGFS.

## Dev subset (Phase 1, `data/interim/dev_subset/`)

- `selection.yaml`: 8 lines × 50 drugs (+DMSO_TF), chosen by
  `scripts/select_dev_subset.py` — 2× KRAS-G12C (HOP62, MIA PaCa-2), 2× KRAS
  non-G12C (SW480 G12V, NCI-H460 Q61H), 2× BRAF-V600E (RKO, LOX-IMVI), 2×
  RAS/RAF-WT (A498, A-172); all 26 fine-MOA classes covered; ~6.08M cells estimated.
- Built by `scripts/build_dev_subset.py` (resumable scan of all shards; adds
  `conc`/`conc_unit` columns; same schema otherwise).
