#!/bin/bash
#SBATCH --job-name=celldrug
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=logs/celldrug_%j.out
cd /mnt/beegfs/scratch/miska/nm667/inProgress/againg/singleCellPerturbationData
export PYTHONPATH=$PWD/src
./.venv/bin/python -m perturbmodel.celldrug --variance
