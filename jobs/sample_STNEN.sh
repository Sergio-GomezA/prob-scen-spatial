# Sample extraction for ST-NEN models
#!/bin/bash
#$ -N STNEN
#$ -wd /exports/eddie/scratch/s2441782/scenarios/prob-scen-spatial/
#$ -o /exports/eddie/scratch/s2441782/scenarios/spatial/ST-NEN/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/scenarios/spatial/ST-NEN/jobfiles/
##$ -l h_rt=4:00:0,h_vmem=16G
#$ -pe sharedmem 8
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -t 1-60

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript model_samples.R $SGE_TASK_ID 1 TRUE ST-NEN r_err.cf_f_gaussian_fd_feat_fcst_group-matern-ar1.rds
