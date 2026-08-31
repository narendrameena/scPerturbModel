# Single-Cell Perturbation Modeling (Tahoe-100M)

Goal: build a predictive model of cellular state based on single-cell perturbation data —
learning how drug perturbations reshape transcriptomes across genetic contexts, using the
**Tahoe-100M** atlas (Zhang et al., bioRxiv 2025, doi:10.1101/2025.02.20.639398) as the
primary training resource. **Framework: PyTorch** (+ Lightning for training loops); scVI
via scvi-tools (itself PyTorch) as the paper-matching baseline.

## Key resources

| Resource | Where |
|---|---|
| Paper (preprint PDF) | `papers/2025.02.20.639398v3.full.pdf` |
| Paper notes / summary | `docs/paper_notes/tahoe-100M.md` |
| Dataset (100M cells) | Hugging Face: `tahoebio/Tahoe-100M` |
| Modeling roadmap | `docs/project_roadmap.md` |
| **Findings so far** | **`RESULTS.md`** |
| Related work survey | `docs/related_work_perturbation_models.md` |
| Methodological positioning | `docs/methodology_positioning.md` |

## Dataset at a glance

- ~100.6M cells passing minimal QC, **95.6M passing full QC** (used for analyses)
- **50 cancer cell lines** (47 analyzed) from 13 organs, mixed as "cell villages" (Mosaic platform)
- **379 distinct drugs** / 1,135 drug–dose combinations → **52,886 unique cell line × drug × dose conditions** (median ~1,287 cells/condition)
- 24 h treatment, Parse GigaLab combinatorial barcoding, SNP-based demultiplexing (demuxlet)
- 14 × 96-well plates; plate 14 is a biological replicate of plate 6 (held out for validation in the paper)
- DMSO vehicle controls on every plate (wells H11/H12) → use for plate-aware differential comparisons

## Directory layout

```
├── papers/            # literature PDFs
├── docs/              # paper notes, roadmap, design docs
├── data/
│   ├── raw/           # untouched downloads (Tahoe-100M shards from HF)
│   ├── external/      # other datasets: sciPlex3, Replogle Perturb-seq, DepMap, scPerturb
│   ├── interim/       # intermediate transforms (filtered/subsampled AnnData, pseudobulk)
│   ├── processed/     # model-ready tensors / train-val-test splits
│   └── metadata/      # drug, MOA, target, cell line, mutation annotation tables
├── notebooks/         # exploratory Jupyter notebooks (numbered: 01_..., 02_...)
├── src/perturbmodel/  # importable Python package (data, preprocessing, models, training, evaluation, utils)
├── scripts/           # CLI entry points (download, preprocess, train, evaluate)
├── configs/           # YAML experiment configs
├── envs/              # conda/pip environment specs
├── jobs/              # SLURM / cluster submission scripts
├── results/           # figures, tables, checkpoints, benchmarks
├── logs/              # run logs (SLURM stdout/stderr, training logs)
└── tests/             # unit tests for src/
```

## Quickstart

```bash
# 1. Create environment
conda env create -f envs/environment.yml
conda activate tahoe

# 2. Inspect what's in the HF dataset repo & fetch metadata tables (small)
python scripts/download_tahoe100m.py --mode list
python scripts/download_tahoe100m.py --mode metadata

# 3. Stream a preview of expression records (no full download needed)
python scripts/download_tahoe100m.py --mode preview -n 5

# 4. (Later, ~hundreds of GB) full snapshot into data/raw/
python scripts/download_tahoe100m.py --mode full
```

## Notes

- This directory lives on BeeGFS **scratch** — fine for the large data, but back up
  code/docs elsewhere (e.g. `git init` + push to a remote).
- Evaluation conventions from the paper worth mirroring: E-distance in a 10-d scVI latent
  space, plate-matched DMSO controls, ≥50 cells per condition, plate 6 vs 14 replication check.
