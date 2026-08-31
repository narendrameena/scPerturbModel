#!/bin/bash
#SBATCH --job-name=cdi_op3
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/drug_cdi_external.py --h5ad data/external/op3/pseudobulk.h5ad \
  --tag OP3 --context cell_type --perturbation sm_name --replicate donor_id \
  --dose dose_uM --control-value "Dimethyl Sulfoxide"
