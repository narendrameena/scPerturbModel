#!/bin/bash
#SBATCH --job-name=dev47_delta
#SBATCH --partition=2204
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/train_delta_model.py --pb-dir data/processed/pseudobulk_dev47 --tag dev47
