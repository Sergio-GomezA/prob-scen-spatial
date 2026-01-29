# Model fit
mstart_t <- Sys.time()

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE

# Get task ID and others from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# browser()
# Assign default values
model_id <- 4
model_list_file <- "data/model_list_spatial.parquet"
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

require(parallel)

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
  mc <- detectCores()
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  .libPaths(temp_lib)
}

require(tidyverse)
require(data.table)
require(lubridate)
require(INLA)
require(arrow)
require(sf)
require(fmesher)
require(ggplot2)
require(ggthemes)

input_data <- "data/scottish_wfsamp_24.parquet"
# input_scalingpars <- "data/scottish_wf_24_scaling_pars.csv"

inla.setOption(num.threads = sprintf("%d:1", mc))
source("fcst_functions.R")
source("functions_probscen.R")
source("aux_funct_ps.R")

t1 <- "2024-06-30 23:00:00" %>% as.POSIXct(tz = "UTC")
# training window in months
window <- 1

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
data_masked <- history_window(
  data.scaled,
  t1,
  window = window,
  mask = mask_opt
) %>%
  mutate(
    site_id = as.integer(factor(site_id)),
    # time_num = as.numeric(time),
    # t = ave(
    #   time_num,
    #   site_id,
    #   FUN = function(x) seq_along(x)
    # )
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_transform(crs = 27700) %>%
  mutate(
    x = st_coordinates(.)[, 1] / 1000,
    y = st_coordinates(.)[, 2] / 1000
  ) %>%
  st_drop_geometry()

time_grid <- seq(
  from = min(data_masked$time),
  to = max(data_masked$time),
  by = "1 hour"
)
data_masked$time_idx <- match(data_masked$time, time_grid)
data_masked$time_idx %>% range()
# A1 <- inla.spde.make.A(
#   mesh = wf.mesh,
#   loc = cbind(data_masked$x, data_masked$y),
#   group = data_masked$t
# )

# wf.spde <- inla.spde2.pcmatern(
#   mesh = wf.mesh,
#   alpha = 2,
#   prior.range = c(20, 0.5), # P(range < 50km) = 0.5
#   prior.sigma = c(1, 0.5)
# )
# spde_idx <- inla.spde.make.index(
#   name = "spatial",
#   n.spde = wf.spde$n.spde,
#   n.group = length(time_grid)
# )
# ncol(A1)
# wf.spde$n.spde * length(time_grid)
# length(spde_idx$spatial)
# data_masked %>%
#   select(x, y) %>%
#   unique() %>%
#   dist()
loc_unique <- data_masked %>%
  distinct(x, y) %>%
  as.matrix()
bnd <- fm_extensions(loc_unique, convex = c(-.2, -.3))
# ggplot() + geom_sf(data = bnd[[1]])
wf.mesh <- fm_mesh_2d(
  # loc = loc_unique,
  loc = fm_hexagon_lattice(bnd[[1]], edge_len = 15),
  boundary = bnd,
  max.edge = c(20, 40), # km
  # offset = -0.2,
  cutoff = 5
)
wf.mesh$n


ggplot() +
  geom_fm(data = wf.mesh) +
  geom_point(
    aes(x, y),
    data = data_masked %>% select(x, y) %>% unique(),
    inherit.aes = FALSE
  ) +
  theme_map()
ggsave("fig/meshhex_scottish_wfsamp.pdf", width = 6, height = 4)
########## Model specifications ###############################################

# read model list
model_list <- read_parquet(model_list_file)
cat(
  sprintf(
    "Model list %s loaded\n",
    model_list_file
  )
)


# extract current choice
model_type <- model_list[model_id, 1:6]
features_vec <- model_list[model_id, 7] %>%
  unlist() %>%
  unname() #%>%
# gsub("ar2", "ar1g", .)
# gsub("ar2", "matern-ar1", .)
# features_vec <- features_vec[-3]
cat(
  sprintf(
    "Running model type: %s \nFeatures included: %s",
    paste0(model_type, collapse = ", "),
    paste0(features_vec, collapse = ", ")
  )
)


########## Model fitting ######################################################
source("aux_funct_ps.R")
# undebug(history_window)
# debug(fit_inla_model)
# initial_values <- NULL

mod_temp <- tryCatch(
  fit_inla_model(
    model_type = model_type,
    features_vec = features_vec,
    data = data_masked,
    ini.theta = initial_values,
    mesh = wf.mesh,
    verbose = TRUE
  ),
  error = function(e) {
    message("Error in model fitting: ", e$message)
    return(NULL)
  }
)

########## Saving output ######################################################

if (!is.null(mod_temp)) {
  # print message
  cat(
    "Model successfully fitted, extracting data ...\n"
  )
  print(summary(mod_temp))
  # extracting hyperparameters' posterior mode
  model_mode <- mod_temp$mode$theta %>% as.list() %>% as.data.frame()
  # extracting hyperparameters' tags
  mode_tags <- mod_temp$mode$theta.tags %>% paste0(., collapse = ", ")
  # extracting scores
  model_scores <- extract_score_model(mod_temp)

  # save_model_info
  result <- bind_cols(model_scores, model_mode, tags = mode_tags)

  model_fname <- sprintf(
    "r_%s_f_%s_%s_feat_%s.rds",
    model_type$response,
    model_type$family,
    model_type$fderiv,
    # model list
    # sub("^model_(.*)\\.rds$", "\\1", model_list_file),
    # id
    paste(features_vec, collapse = "-")
  )
  if (save_model) saveRDS(mod_temp, file = file.path(mod_obj_path, model_fname))
} else {
  cat("Model fitting step failed \n")
  result <- data.frame()
}

# store result in output folder
# create directories if necessary
if (!dir.exists(path_out)) {
  dir.create(path_out, recursive = TRUE)
}

fname <- file.path(
  path_out,
  sprintf(
    "var_scores_r_%s_f_%s_%s_feat_%s.csv",
    model_type$response,
    model_type$family,
    model_type$fderiv,
    # model list
    # sub("^model_(.*)\\.rds$", "\\1", model_list_file),
    # id
    paste(features_vec, collapse = "-")
  )
)
write.csv(result, fname, row.names = FALSE)
# print message
cat(sprintf("Saving output in file %s\n", fname))

mend_t <- Sys.time()
run_time <- difftime(mend_t, mstart_t)
cat(sprintf(
  "Whole process ended for model id: %d in %.2f %s\n",
  model_id,
  run_time,
  units(run_time)
))
