#!/bin/bash
#SBATCH --job-name=full_var
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/decompose_response_variance.py --pb-dir data/processed/pseudobulk_full --tag full
./.venv/bin/python scripts/diagnose_context_structure.py --pb-dir data/processed/pseudobulk_full --tag full
