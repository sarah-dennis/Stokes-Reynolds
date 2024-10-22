#!/bin/bash
#SBATCH --account=guest
#SBATCH --partition=guest-compute
#SBATCH --job-name=hex-Re05
#SBATCH --ntasks=1
#SBATCH --output=_log_%j.txt

module load share_modules/ANACONDA/5.3_py3
python3 stokes_control_run.py $SLURM_ARRAY_TASK_ID