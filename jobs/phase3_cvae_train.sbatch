#!/bin/bash
#SBATCH --job-name=cvae_train
#SBATCH --partition=2204
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
echo "Job ${SLURM_JOB_ID} on $(hostname), started $(date)"
./.venv/bin/python scripts/train_cond_vae.py
echo "Finished $(date)"
