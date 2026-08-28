#!/bin/bash
#SBATCH --job-name=score_cells
#SBATCH --partition=2204
#SBATCH --cpus-per-task=20
#SBATCH --mem=300G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/score_cells_on_components.py --src data/interim/dev47/data --tag dev47 --workers 16
