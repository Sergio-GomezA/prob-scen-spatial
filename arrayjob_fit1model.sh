#!/bin/bash

#$ -N psST
#$ -wd /exports/eddie/scratch/s2441782/scenarios/prob-scen-spatial/
#$ -o /exports/eddie/scratch/s2441782/scenarios/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/scenarios/jobfiles/
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -pe sharedmem 8
#$ -t 1-4

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run model fitting code
Rscript model_fit.R $SGE_TASK_ID
