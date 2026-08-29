#!/bin/bash
#SBATCH --job-name=drug_ctx
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/drug_context_dependence.py --pb-dir data/processed/pseudobulk_full --tag full
