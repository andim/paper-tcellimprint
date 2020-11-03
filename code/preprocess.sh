#!/bin/bash 

# set number of nodes
#SBATCH -N 1
#SBATCH --ntasks-per-node=1 
#SBATCH --mem=8000MB 
# time in format dd-hh:mm:ss
#SBATCH --time=01-00:00:00
# set nice descriptive name 
#SBATCH -J preprocess
# queue name
#SBATCH -p bioth
#SBATCH --array=1-786
#SBATCH --output=test.out

source activate py3
python -V

echo $SLURM_ARRAY_TASK_ID $HOSTNAME 
python preprocess_emerson.py $SLURM_ARRAY_TASK_ID
