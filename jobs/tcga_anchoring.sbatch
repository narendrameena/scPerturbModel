#!/bin/bash
#SBATCH --job-name=tcga
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/tcga_anchoring.py --tag dev47
