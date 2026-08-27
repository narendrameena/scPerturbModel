#!/bin/bash
#SBATCH --job-name=ctx_struct
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/diagnose_context_structure.py --pb-dir data/processed/pseudobulk_dev47 --tag dev47
