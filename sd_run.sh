#!/bin/bash
#SBATCH --account=guest
#SBATCH --partition=guest-compute
#SBATCH --job-name=hex-Re0-Re1
#SBATCH --output=_log_%j.txt
#SBATCH --ntasks=1

module load share_modules/ANACONDA/5.3_py3
python3 stokes_control_run.py $SLURM_ARRAY_TASK_ID