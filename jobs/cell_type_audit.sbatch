#!/bin/bash
#SBATCH --job-name=ct_audit
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
echo "Job ${SLURM_JOB_ID} on $(hostname), started $(date)"
./.venv/bin/python scripts/audit_cell_types.py
echo "Finished $(date)"
