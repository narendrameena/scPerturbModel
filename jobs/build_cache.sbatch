#!/bin/bash
#SBATCH --job-name=cache
#SBATCH --partition=2204
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=logs/cache_%j.out
cd /mnt/beegfs/scratch/miska/nm667/inProgress/againg/singleCellPerturbationData
export PYTHONPATH=$PWD/src
./.venv/bin/python -c "
from perturbmodel.celldrug import prism_gamma, apportion
g, ti, dec = prism_gamma(rebuild=True)
print(f'cached {len(g)} compounds')
v = apportion(dec)
import json; print(json.dumps({k: v[k] for k in v if k != 'gamma_covs'}, indent=1))
"
