# Preparing initial dataset of scottish wind farms.

require(tidyverse)
require(data.table)
require(parallel)

require(arrow)

require(rnaturalearth)
require(rnaturalearthdata)
require(sf)

require(ggthemes)
require(ggsci)


# Global settings ####

# read Data ####
if (grepl("exports", getwd())) {
  # running in cluster
  cluster_run <- TRUE
  data_path <- gen_path <- "data"
  model_path <- "../calibration/model_objects"
  gen_path <- "../calibration/data"
  temp_lib <- "/exports/eddie/scratch/s2441782/calibration/lib"
  .libPaths(temp_lib)
  setwd("/exports/eddie/scratch/s2441782/calibration_power")
} else {
  cluster_run <- FALSE
  data_path <- "~/Documents/ERA5_at_wf/"
  gen_path <- "~/Documents/elexon/"
  model_path <- "~/Documents/elexon/model_objects"
}


if (cluster_run) {
  n.cores <- detectCores()
  pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve.parquet"))
} else {
  n.cores <- detectCores() - 2
  pwr_curv_df <- read_parquet(file.path(gen_path, "power_curve.parquet"))
}

source("aux_funct.R")
# inla.setOption(num.threads = paste0(n.cores, ":1"))
