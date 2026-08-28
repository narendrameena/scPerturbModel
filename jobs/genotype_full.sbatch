#!/bin/bash
#SBATCH --job-name=geno_full
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
echo "=== MOA-aggregated scan (primary) ==="
./.venv/bin/python scripts/genotype_drug_interaction.py --pb-dir data/processed/pseudobulk_full --tag fullmoa --by-moa
echo "=== per-drug scan (secondary, 379 drugs) ==="
./.venv/bin/python scripts/genotype_drug_interaction.py --pb-dir data/processed/pseudobulk_full --tag fulldrug
