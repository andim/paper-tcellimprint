#!/bin/bash 

# set number of nodes
#SBATCH -N 1
#SBATCH --ntasks-per-node=1 
#SBATCH --mem=8000MB 
# time in format dd-hh:mm:ss
#SBATCH --time=24-00:00:00
# set nice descriptive name 
#SBATCH -J mayer
# queue name
#SBATCH -p bioth
# run as an array job (change number of tasks here)
#SBATCH --array=1-3
#SBATCH --output=test.out

echo $SLURM_ARRAY_TASK_ID $HOSTNAME 
python run-death.py $SLURM_ARRAY_TASK_ID
