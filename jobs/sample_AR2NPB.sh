# Sample extraction for AR2-NPB models
#!/bin/bash
#$ -N AR2NPB
#$ -wd /exports/eddie/scratch/s2441782/scenarios/prob-scen-spatial/
#$ -o /exports/eddie/scratch/s2441782/scenarios/spatial/AR2-NPB/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/scenarios/spatial/AR2-NPB/jobfiles/
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
Rscript model_samples.R $SGE_TASK_ID 1 TRUE AR2-NPB r_actuals.cf_f_beta_fd_feat_ws.w_group-fcst_group-fd_group-hour-ar2g.rds
