#!/bin/bash
#SBATCH --job-name=probe_panel
#SBATCH --partition=2204
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/eval_few_shot.py --pb-dir data/processed/pseudobulk_dev47 \
  --tag panel --adapt finetune --seeds 1 --k-values 5 20 \
  --strategies random high_effect discriminative moa_diverse chem_diverse
