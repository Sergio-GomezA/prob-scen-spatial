#!/bin/bash

#$ -N ps-st
#$ -wd /exports/eddie/scratch/s2441782/spatial/
#$ -o /exports/eddie/scratch/s2441782/spatial/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/spatial/jobfiles/
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -pe sharedmem 64

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript model_fit.R
