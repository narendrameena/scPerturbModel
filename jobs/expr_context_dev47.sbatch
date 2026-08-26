#!/bin/bash
#SBATCH --job-name=exprctx_dev47
#SBATCH --partition=2204
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/eval_hard_splits.py --mode line --context expression --pb-dir data/processed/pseudobulk_dev47 --tag dev47expr
