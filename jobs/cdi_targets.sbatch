#!/bin/bash
#SBATCH --job-name=cdi_tgt
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/cdi_vs_target_genetics.py --pb-dir data/processed/pseudobulk_full
