# Sample extraction for ST-PR-NPB models
#!/bin/bash
#$ -N STPRNPB
#$ -wd /exports/eddie/scratch/s2441782/scenarios/prob-scen-spatial/
#$ -o /exports/eddie/scratch/s2441782/scenarios/spatial/ST-PR-NPB/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/scenarios/spatial/ST-PR-NPB/jobfiles/
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
Rscript model_samples.R $SGE_TASK_ID 1 TRUE ST-PR-NPB r_actuals.cf_f_beta_eta_feat_fcst_group-matern-ar1-etaderiv.rds
