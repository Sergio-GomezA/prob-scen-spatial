# Sample extraction for AR2-PR-NEN models
#!/bin/bash
#$ -N AR2PRNEN
#$ -wd /exports/eddie/scratch/s2441782/scenarios/prob-scen-spatial/
#$ -o /exports/eddie/scratch/s2441782/scenarios/spatial/AR2-PR-NEN/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/scenarios/spatial/AR2-PR-NEN/jobfiles/
##$ -l h_rt=4:00:0,h_vmem=16G
#$ -pe sharedmem 8
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -t 1-2

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript model_samples.R $SGE_TASK_ID 1 TRUE AR2-PR-NEN r_err.cf_f_gaussian_eta_feat_ws.w_group-fcst_group-hour-ar2g-etaderiv.rds
