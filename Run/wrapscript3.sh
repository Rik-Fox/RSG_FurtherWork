#!/bin/bash
#SBATCH --job-name=HAT_RSG
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=5g
#SBATCH --export=ALL
#SBATCH --array=1-5

mkdir -p tmp/$SLURM_JOB_ID

module --ignore-cache load "matlab"

matlab -nodisplay -nosplash -r "Wrapper3($SLURM_ARRAY_TASK_ID); exit"

rm -rf tmp/$SLURM_JOB_ID
