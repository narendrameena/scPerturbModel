#!/bin/bash
#SBATCH --job-name=epigenome
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/epigenome_vs_context.py --pb-dir data/processed/pseudobulk_full
