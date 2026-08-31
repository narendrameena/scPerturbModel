#!/bin/bash
#SBATCH --job-name=lincs_p1
#SBATCH --partition=2204
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
./.venv/bin/python scripts/lincs_mechanism_ranking.py \
  --gctx data/external/lincs/p1_level4.gctx \
  --gz data/external/lincs/p1_level4.gctx.gz \
  --inst-info data/external/lincs/p1_inst_info.txt.gz \
  --gene-info data/external/lincs/p1_gene_info.txt.gz \
  --tag phase1 --min-lines 10
