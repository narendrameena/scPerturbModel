#!/bin/bash
#SBATCH --job-name=fewshot_dev8
#SBATCH --partition=2204
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/eval_few_shot.py --pb-dir data/processed/pseudobulk_dev --tag dev8 --seeds 3
