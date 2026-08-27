#!/bin/bash
#SBATCH --job-name=fewshot_pilot
#SBATCH --partition=2204
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/eval_few_shot.py --folds 2 --seeds 2 --tag pilot
