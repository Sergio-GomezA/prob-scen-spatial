# Sample extraction for ST-PR-NEN models
#!/bin/bash
#$ -N STPRNEN
#$ -wd /exports/eddie/scratch/s2441782/scenarios/prob-scen-spatial/
#$ -o /exports/eddie/scratch/s2441782/scenarios/spatial/ST-PR-NEN/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/scenarios/spatial/ST-PR-NEN/jobfiles/
##$ -l h_rt=4:00:0,h_vmem=16G
#$ -pe sharedmem 8
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -t 4-60

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript model_samples.R $SGE_TASK_ID 1 TRUE ST-PR-NEN r_err.cf_f_gaussian_eta_feat_ws.w_group-matern-ar1-etaderiv.rds
