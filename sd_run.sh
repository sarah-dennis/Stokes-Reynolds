#!/bin/bash
#SBATCH --account=guest
#SBATCH --partition=guest-compute
#SBATCH --job-name=BFS_Re0
#SBATCH --ntasks=1
#SBATCH --output=_log_%j.txt

module load share_modules/ANACONDA/5.3_py3
python3 stokes_run.py $SLURM_ARRAY_TASK_ID
