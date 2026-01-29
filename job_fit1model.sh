#!/bin/bash

#$ -N psst1
#$ -wd /exports/eddie/scratch/s2441782/spatial/
#$ -o /exports/eddie/scratch/s2441782/spatial/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/spatial/jobfiles/
#$ -M s2441782@ed.ac.uk
#$ -m bea

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run model fitting code
Rscript model_fit.R
