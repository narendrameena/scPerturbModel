#!/bin/bash
#SBATCH --job-name=fewshot_dev47
#SBATCH --partition=2204
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/eval_few_shot.py --seeds 3 --tag dev47
