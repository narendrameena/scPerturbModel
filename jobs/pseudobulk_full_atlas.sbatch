#!/bin/bash
#SBATCH --job-name=pb_full
#SBATCH --partition=2204
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/build_pseudobulk.py --raw \
  --src data/raw/tahoe-100m/data --out data/processed/pseudobulk_full
