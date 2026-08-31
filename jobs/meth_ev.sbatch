#!/bin/bash
#SBATCH --job-name=meth_ev
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=10:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/methodology_evidence.py
