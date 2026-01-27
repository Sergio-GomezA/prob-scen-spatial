# Model fit

require(tidyverse)
require(data.table)
require(lubridate)
require(INLA)
require(arrow)

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

# Get task ID and others from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# browser()
# Assign default values
model_id <- 4
model_list_file <- "data/model_top_eta2.rds"
ofolder <- "etaderiv"
save_model <- TRUE
initial_values_file <- ""

# Override defaults only if arguments are provided
if (length(args) > 0) {
  model_id <- as.numeric(args[1])
}
if (length(args) > 1) {
  model_list_file <- as.character(args[2])
}
if (length(args) > 2) {
  ofolder <- as.character(args[3])
}
if (length(args) > 3) {
  save_model <- as.logical(args[4])
}
if (length(args) > 4) {
  initial_values_file <- as.character(args[5])
}

if (local_run) {
  large_obj_path <- "~/Documents/proj2/spatial"
  mod_obj_path <- "~/Documents/proj2/spatial/model_objects"
  path_out <- file.path("~/Documents/proj2", "spatial", ofolder)
  # path.fig <- file.path("~/Documents/proj2",ofolder,"fig")
  # path.samples <- file.path("~/Documents/proj2",ofolder,"sample")
  mc <- detectCores() - 2
} else {
  large_obj_path <- "/exports/eddie/scratch/s2441782/scenarios/spatial"
  mod_obj_path <- file.path(large_obj_path, "model_objects")
  path_out <- file.path(
    "/exports/eddie/scratch/s2441782/scenarios/spatial",
    ofolder
  )
  # path.fig <- file.path(large_obj_path,ofolder,"fig")
  # path.samples <- file.path(large_obj_path,ofolder,"sample")
  mc <- 1 #detectCores()-2
}

input_data <- "data/scotish_wf_24.parquet"
# input_scalingpars <- "data/scotish_wf_24_scaling_pars.csv"

inla.setOption(num.threads = sprintf("%d:1", mc))
source("fcst_functions.R")
source("functions_probscen.R")
source("aux_funct_ps.R")

t1 <- "2024-01-01 00:00:00" %>% as.POSIXct(tz = "UTC")
# training window in months
window <- 6

# remove last 24 hours
mask_opt <- TRUE

# initial values for hyper parameters

if (initial_values_file != "") {
  # read top scores
  top_scores <- readRDS(initial_values_file)
  # get initials
  initial_values <- top_scores[[model_id]]
} else {
  initial_values <- list()
}

cens <- 0.001

# priors
hyper.rw2 <- list(
  theta = list(initial = 1, param = c(0.2, 0.0001), fixed = FALSE)
)
hyper.ar1 <- list(
  theta1 = list(initial = 1, param = c(0.2, 0.0001), fixed = FALSE), #logprecision
  theta2 = list(initial = 0, param = c(0, 0.15), fixed = FALSE) # rho
)

hyper.ar2 <- list(
  theta1 = list(initial = 1, param = c(0.2, 0.0001), fixed = FALSE), #log-precision
  theta2 = list(initial = 3, param = c(0.5, 0.5), fixed = FALSE), # pacf1
  theta3 = list(initial = -1, param = c(0.4, 0.4), fixed = FALSE) # pacf1
)
hyper.iid <- list(theta = list(initial = 5, param = c(1, 0.01), fixed = FALSE))

# mycontrol_inla <- control.inla(control.vb = control.vb(enable = FALSE, strategy = "mean"))

if (!dir.exists(mod_obj_path)) {
  dir.create(mod_obj_path, recursive = TRUE)
}

########## Data loading #######################################################
# input_data <- "data/aggr_powr_bpa_12-23.csv.gz"
data.scaled <- read_parquet(input_data)
# scaling_params <- read_parquet(input_scalingpars)
cat(sprintf(
  "Input file %s loaded\n",
  input_data
))


########## Model specifications ###############################################

# read model list
model_list <- readRDS(model_list_file)
cat(
  sprintf(
    "Model list %s loaded\n",
    model_list_file
  )
)

# extract current choice
model_type <- model_list[model_id, 1:5]
features_vec <- model_list[model_id, 6] %>% unlist() %>% unname()

cat(
  sprintf(
    "Running model type: %s \nFeatures included: %s",
    paste0(model_type, collapse = ", "),
    paste0(features_vec, collapse = ", ")
  )
)

########## Model fitting ######################################################

mod_temp <- tryCatch(
  fit_inla_model(
    model_type = model_type,
    features_vec = features_vec,
    data = history_window(
      data.scaled,
      t1,
      window = window,
      mask = mask_opt
    ),
    ini.theta = initial_values,
    verbose = TRUE
  ),
  error = function(e) {
    message("Error in model fitting: ", e$message)
    return(NULL)
  }
)
