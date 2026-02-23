# Sample extraction for AR2-PR-NPB_f models
#!/bin/bash
#$ -N AR2PRNPB_f
#$ -wd /exports/eddie/scratch/s2441782/scenarios/prob-scen-spatial/
#$ -o /exports/eddie/scratch/s2441782/scenarios/spatial/AR2-PR-NPB_f/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/scenarios/spatial/AR2-PR-NPB_f/jobfiles/
##$ -l h_rt=4:00:0,h_vmem=16G
#$ -pe sharedmem 8
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -t 3-60

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript model_samples.R $SGE_TASK_ID 1 TRUE AR2-PR-NPB_f r_actuals.cf_f_beta_eta_feat_fcst_group-ar2g-etaderiv.rds
