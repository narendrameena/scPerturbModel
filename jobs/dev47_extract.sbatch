#!/bin/bash
#SBATCH --job-name=dev47_extract
#SBATCH --partition=2204
#SBATCH --cpus-per-task=28
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/build_dev_subset.py --subset-dir data/interim/dev47 --workers 24
