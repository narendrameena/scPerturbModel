#!/bin/bash
#SBATCH --job-name=geno_drug
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/genotype_drug_interaction.py --pb-dir data/processed/pseudobulk_dev47 --tag dev47
