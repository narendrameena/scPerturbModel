#!/bin/bash
#SBATCH --job-name=tahoe
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
##SBATCH --partition=<PARTITION>          # <- set for this cluster
##SBATCH --gres=gpu:1                     # <- uncomment for GPU jobs
##SBATCH --account=<ACCOUNT>

# Usage: sbatch jobs/template.sbatch scripts/<script>.py --config configs/default.yaml
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

# module load cuda/12.x   # if the cluster uses modules
source activate tahoe 2>/dev/null || conda activate tahoe

echo "Job $SLURM_JOB_ID on $(hostname), started $(date)"
python "$@"
echo "Finished $(date)"
