#!/bin/bash
#SBATCH --job-name=dev47_pb
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=16:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/build_pseudobulk.py --src data/interim/dev47/data --out data/processed/pseudobulk_dev47
