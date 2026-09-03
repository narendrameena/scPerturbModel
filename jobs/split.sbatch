#!/bin/bash
#SBATCH --job-name=split
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=logs/celldrug_%j.out
cd /mnt/beegfs/scratch/miska/nm667/inProgress/againg/singleCellPerturbationData
export PYTHONPATH=$PWD/src
./.venv/bin/python scripts/state_vs_genotype_split.py --n-compounds 150
