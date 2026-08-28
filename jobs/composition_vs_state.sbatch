#!/bin/bash
#SBATCH --job-name=comp_state
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/composition_vs_state.py --tag dev47
