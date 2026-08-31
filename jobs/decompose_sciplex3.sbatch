#!/bin/bash
#SBATCH --job-name=dec_sp3
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/pseudobulk_sciplex3.py
./.venv/bin/python scripts/decompose_external.py \
  --h5ad data/external/scperturb/sciplex3_pseudobulk.h5ad \
  --tag sciPlex3 --context cell_line --perturbation perturbation \
  --replicate replicate --dose dose_value --control-value control
