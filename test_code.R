# pcov <- myspde.posterior(mod.temp, "spatial", "matern.covariance")
# # pcov %>% names()
# pcov %>%
#   ggplot() +
#   geom_ribbon(
#     aes(x, ymin = `q0.025%`, ymax = `q0.975%`),
#     alpha = 0.5
#   ) +
#   gg(pcov)

# pcor <- myspde.posterior(mod.temp, "spatial", "matern.correlation")
# pcor %>%
#   ggplot() +
#   geom_ribbon(
#     aes(x, ymin = `q0.025%`, ymax = `q0.975%`),
#     alpha = 0.5
#   ) +
#   gg(pcor)

model_fnames <- list.files("~/Documents/proj2/spatial/model_objects")

test_df <- data.frame(
  mod.file.name = model_fnames
) %>%
  mutate(
    ofolder = case_when(
      grepl("err.cf", mod.file.name) &
        grepl("etaderiv", mod.file.name) &
        grepl("matern", mod.file.name) ~ "ST-PR-NEN",
      grepl("actuals.cf", mod.file.name) &
        grepl("etaderiv", mod.file.name) &
        grepl("matern", mod.file.name) ~ "ST-PR-NPB",
      grepl("err.cf", mod.file.name) &
        grepl("etaderiv", mod.file.name) &
        grepl("ar1g", mod.file.name) ~ "AR1-PR-NEN",
      grepl("actuals.cf", mod.file.name) &
        grepl("etaderiv", mod.file.name) &
        grepl("ar1g", mod.file.name) ~ "AR1-PR-NPB",
      grepl("err.cf", mod.file.name) &
        grepl("matern", mod.file.name) ~ "ST-NEN",
      grepl("actuals.cf", mod.file.name) &
        grepl("matern", mod.file.name) ~ "ST-NPB",
      grepl("err.cf", mod.file.name) &
        grepl("etaderiv", mod.file.name) ~ "PR-NEN",
      grepl("actuals.cf", mod.file.name) &
        grepl("etaderiv", mod.file.name) ~ "PR-NPB",
      TRUE ~ "Other"
    )
  )

scores_path <- "~/Documents/proj2/spatial/etaderiv"
cpo_scores <- list.files(
  scores_path
) %>%
  lapply(., \(x) fread(file.path(scores_path, x)) %>% mutate(fname = x)) %>%
  bind_rows() %>%
  mutate(
    mod.file.name = paste0(gsub("var_scores_|\\.csv", "", fname), ".rds")
  ) %>%
  right_join(., test_df, by = "mod.file.name") %>%
  select(ofolder, mod.file.name, everything()) %>%
  group_by(
    ofolder
  ) %>%
  filter(
    !grepl("ws.w_group-fcst_group", mod.file.name)
  ) %>%
  filter(
    waic == min(waic)
  )

write.csv(cpo_scores, "summaries/top_stmodel_cpo.csv", row.names = F)


##

# Obtaining samples from models

# job files

"# Sample extraction for ST-PR-NEN models
#!/bin/bash
#$ -N STPRNEN
#$ -wd /exports/eddie/scratch/s2441782/scenarios/prob-scen-spatial/
#$ -o /exports/eddie/scratch/s2441782/scenarios/spatial/ST-PR-NEN/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/scenarios/spatial/ST-PR-NEN/jobfiles/
##$ -l h_rt=4:00:0,h_vmem=16G
#$ -pe sharedmem 8
#$ -M s2441782@ed.ac.uk
#$ -m bea
#$ -t 1-3

# Initialise modules
source /etc/profile.d/modules.sh

# Load R
module load R/4.5

# Run resolution code
Rscript model_samples.R $SGE_TASK_ID 1 TRUE ST-PR-NEN r_err.cf_f_gaussian_eta_feat_ws.w_group-matern-ar1-etaderiv.rds
"

#
for (i in seq_along(cpo_scores$mod.file.name)) {
  # browser()
  jobname <- cpo_scores$ofolder[i] %>% gsub("-", "", .)

  job_script <- sprintf("jobs/sample_%s.sh", jobname)

  script_content <- sprintf(
    "# Sample extraction for %s models
#!/bin/bash
#$ -N %s
#$ -wd /exports/eddie/scratch/s2441782/scenarios/prob-scen-spatial/
#$ -o /exports/eddie/scratch/s2441782/scenarios/spatial/%s/jobfiles/
#$ -e /exports/eddie/scratch/s2441782/scenarios/spatial/%s/jobfiles/
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
Rscript model_samples.R $SGE_TASK_ID 1 TRUE %s %s",
    cpo_scores$ofolder[i], # dash name description
    jobname, #  no-dash name
    cpo_scores$ofolder[i], # ofolder
    cpo_scores$ofolder[i], # ofolder
    cpo_scores$ofolder[i], # ofolder
    cpo_scores$mod.file.name[i] # file name
  )
  # Write the script to a file
  writeLines(script_content, con = job_script)

  # Make the script executable
  Sys.chmod(job_script, mode = "0755")
}
