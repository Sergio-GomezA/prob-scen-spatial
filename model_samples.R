# author: Sergio

# Probabilistic scenarios script for running in CLI

mstart_t <- Sys.time()
########## Loading packages ###################################################

require(tidyverse)
require(parallel)
require(sf)

########## Global variables ###################################################

local_run <- if (startsWith(getwd(), "/home/s2441782")) TRUE else FALSE
# local_run <- FALSE

# job id is day index to get samples
# Defaults if no argument is provided
day_id <- 25
window_id <- 1
restart <- TRUE
# mod.file.name <- "2503_us_npow_eta0.rds"
mod.file.name <- "r_err.cf_f_gaussian_eta_feat_fcst_group-matern-ar1-etaderiv.rds"
mod.file.name <- "r_actuals.cf_f_beta_eta_feat_fcst_group-matern-ar1-etaderiv.rds"
mod.file.name <- "r_err.cf_f_gaussian_eta_feat_fcst_group-ar1g-etaderiv.rds"
ofolder <- case_when(
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

# Get task ID and others from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
regime <- NULL #"mid"

if (length(args) > 0) {
  day_id <- as.numeric(args[1])
}
if (length(args) > 1) {
  window_id <- as.numeric(args[2])
} # Convert argument to number
if (length(args) > 2) {
  restart <- as.logical(args[3])
}
if (length(args) > 3) {
  ofolder <- as.character(args[4])
}
if (length(args) > 4) {
  mod.file.name <- as.character(args[5])
}
if (length(args) > 5) {
  regime <- as.character(args[6])
}

mod_suffix <- str_extract(mod.file.name, "(?<=feat_).*?(?=\\.rds$)")
# mod.version <- str_extract(mod.file.name, "[^r_.]+(?=\\.rds$)")
mod.version <- str_extract(mod.file.name, "(?<=reg_).*?(?=\\.rds$)")

hregime <- case_when(
  mod.version %in% c("l") ~ "low20%",
  mod.version %in% c("m") ~ "mid",
  mod.version %in% c("h") ~ "high20%",
  TRUE ~ "NULL"
)
hregime <- if (hregime == "NULL") NULL else hregime

if (local_run) {
  large_obj_path <- "~/Documents/proj2/spatial"
  # input_scalingpars <- "aggr_powr_bpa_scaling_pars.csv"
} else {
  large_obj_path <- "/exports/eddie/scratch/s2441782/scenarios/spatial"
  temp_lib <- "/exports/eddie3_homes_local/s2441782/lib"
  .libPaths(temp_lib)
}

require(ggsci)

# require(data.table)
require(arrow)
require(lubridate)
require(stringr)
require(INLA)
require(fmesher)

input_data <- "data/scottish_wfsamp_24.parquet"
path.fig <- file.path(large_obj_path, ofolder, "fig")
path.samples <- file.path(large_obj_path, ofolder, "sample")
mod_obj_path <- file.path(
  large_obj_path,
  ifelse(is.null(hregime), "model_objects", "model_objects_reg")
)

source("fcst_functions.R")
source("functions_probscen.R")
source("aux_funct.R")
source("aux_funct_ps.R")

mc <- available_cores() - ifelse(local_run, 2, 0)
inla.setOption(num.threads = sprintf("%d:1", mc))

# training data window length in units
window_lengths <- c(7, 14, 21, 30, 2)
window_units <- c(rep("days", 4), "months")

# first day to forecast
t0 <- as.POSIXct("2024-07-01", tz = "UTC")

# days to run
time_seq <- seq(
  from = t0,
  to = t0 + 60 * 24 * 3600,
  by = "day"
)

# hours to predict
h_predict <- 24

quantile.seq <- c(0.025, 0.5, 0.975)
n.samples <- 1000
show.fig <- FALSE
save.fig <- TRUE
cens <- 0.001

# translate variables from arguments
t1 <- time_seq[day_id]
window <- window_lengths[window_id]
units <- window_units[window_id]

# Priors
hyper.rw2 <- list(theta = list(initial = 5, param = c(1, 0.01), fixed = FALSE))
hyper.ar1 <- list(
  theta1 = list(initial = 5, param = c(1, 5e-05), fixed = FALSE), #log-precision
  theta2 = list(initial = 3, param = c(0, 0.15), fixed = FALSE) # rho
)
hyper.ar2 <- list(
  theta1 = list(initial = 5, param = c(3, 0.01), fixed = FALSE), #log-precision
  theta2 = list(initial = 3, param = c(0.5, 0.5), fixed = FALSE), # pacf1
  theta3 = list(initial = -1, param = c(0.4, 0.4), fixed = FALSE) # pacf1
)
hyper.iid <- list(theta = list(initial = 5, param = c(1, 0.01), fixed = FALSE))

resp_columns <- c("actuals.cf", "actuals", "err.cf", "err")
group_var <- "power.group"

########## Data loading #######################################################
# input_data <- "data/aggr_powr_bpa_12-23.csv.gz"
data.scaled <- read_parquet(input_data)
cat(sprintf(
  "
Input file %s loaded\nNumber of rows in full data: %d\nNumber of wind farms: %d
Spanning dates: %s to %s
",
  input_data,
  nrow(data.scaled),
  data.scaled %>% pull(site_name) %>% unique() %>% length(),
  min(data.scaled$time) %>% format(., "%Y-%m-%d %H:%M"),
  max(data.scaled$time) %>% format(., "%Y-%m-%d %H:%M")
))

# regime groups consitent list of dates
if (!is.null(regime)) {
  # t0 <- min(time_seq)
  regime_names <- c("low20%", "mid", "high20%")
  day_filters <-
    lapply(
      regime_names,
      \(regime) {
        with(
          data.scaled %>%
            group_by(date) %>%
            summarise(power.group = first(power.group)) %>%
            filter(date >= t0, date <= t0 + 183 * 24 * 3600),
          which(power.group %in% c(regime))
        )
      }
    ) %>%
    setNames(regime_names)
  t1 <- time_seq[day_filters[[regime]]][day_id]
}

########## Filtering data #####################################################

train_data <- history_window(
  data.scaled,
  t1,
  window = window,
  units = units,
  mask = TRUE
) %>%
  mutate(
    site_id = as.integer(factor(site_id))
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
  st_drop_geometry() %>%
  # mask response in other regimes
  {
    if (!is.null(hregime)) {
      mutate(
        .,
        across(
          any_of(resp_columns),
          \(x) ifelse(.data[[group_var]] %in% hregime, x, NA)
        )
      )
    } else {
      .
    }
  }

time_grid <- seq(
  from = min(train_data$time),
  to = max(train_data$time),
  by = "1 hour"
)

train_data$time_idx <- match(train_data$time, time_grid)
time_idx_range <- train_data$time_idx %>% range()
cat(sprintf(
  "Time reindexed ranging from %d to %d hours after initial time\n",
  time_idx_range[1],
  time_idx_range[2]
))

cat("Building spatial mesh\n")
loc_unique <- train_data %>%
  distinct(x, y) %>%
  as.matrix()
bnd <- fm_extensions(loc_unique, convex = c(-.15, -.25))
# bnd <- fm_extensions(loc_unique, convex = c(-.1, -.15))
# ggplot() + geom_sf(data = bnd[[1]])
wf.mesh <- fm_mesh_2d(
  # loc = loc_unique,
  loc = fm_hexagon_lattice(bnd[[1]], edge_len = 15),
  boundary = bnd,
  max.edge = c(20, 40), # km
  # offset = -0.2,
  cutoff = 5
)

# wf.mesh$n

# ggplot() +
#   geom_fm(data = wf.mesh) +
#   geom_point(
#     aes(x, y),
#     data = data_masked %>% select(x, y) %>% unique(),
#     inherit.aes = FALSE
#   ) +
#   theme_map()

cat(sprintf(
  "Data loaded and filtered to %d records using options
Initial date: %s
Window: %d %s\n",
  nrow(train_data),
  format(t1, "%Y-%m-%d"),
  window,
  ifelse(window == 1, sub("s", "", units), units)
))


########## Model fit ##########################################################
cat(sprintf('Reading model file\n'))
mod.path <- file.path(mod_obj_path, mod.file.name)
mod.temp <- readRDS(mod.path)
mod_stats <- get_mod_stats(mod.file.name, path = mod_obj_path, mod.temp)

cat(
  sprintf(
    "Model object %s read\n",
    mod.file.name
  )
)

response <- mod_stats$response.code
modfamily <- mod_stats$family

mod.code <- sprintf(
  "r_%s_f_%s_%s_feat_%s_t",
  response %>% sub(pattern = "Y_", "", .),
  modfamily,
  ifelse(mod_stats$etaderiv, "eta", "fd"),
  mod_suffix
) #%>%
# gsub("#\\.", "-", .)

# Refit model with filtered data
cat("Refitting model for new day-ahead prediction")
# record initial time
start_time <- Sys.time()
source("aux_funct_ps.R")
# debug(fit_a_date)
# re fit model
new_fit <- fit_a_date(
  timet = t1,
  h = h_predict,
  dat.power = train_data,
  response = response,
  inla.object = mod.temp,
  restart = restart,
  mesh = wf.mesh,
  save_stack = TRUE,
  verbose = TRUE
)
summary(new_fit)
# record initial time
end_time <- Sys.time()
run_time <- difftime(end_time, start_time)
cat(sprintf("Refit step finished in: %.2f %s\n", run_time, units(run_time)))


########## Extract samples ####################################################
cat('Extracting posterios samples from model\n')
# record initial time
start_time <- Sys.time()
# labels according to model response
resp.lab <- case_when(
  grepl("actuals", response) ~ "Wind Generation",
  grepl("err", response) ~ "Forecast error"
)

# create directories if necessary
if (!dir.exists(path.samples)) {
  dir.create(path.samples, recursive = TRUE)
}
if (!dir.exists(path.fig)) {
  dir.create(path.fig, recursive = TRUE)
}
# undebug(simulation.plots.inla2)
# save samples plot and get posterior samples

set.seed(1)
source("aux_funct_ps.R")
# debug(simulation.plots.inla2)
# debug(plot_actuals_model)
sim.obj <- simulation.plots.inla2(
  inla.model = new_fit,
  data = data.scaled,
  response = response,
  t1 = t1,
  quSeq = quantile.seq,
  family = modfamily,
  resp.lab = resp.lab,
  # ylim = c(0,1.02),
  nsamp = n.samples,
  show.fig = show.fig,
  save.fig = save.fig,
  fig.ext = ".png",
  path = path.fig,
  run.name = paste0(mod.code, t1),
  skip.plots = FALSE,
  inla_seed = 1
  # sample.df = sample.test$samples,
  # legend.position = "bottom"
)

end_time <- Sys.time()
run_time <- difftime(end_time, start_time)
cat(sprintf(
  "Sample extraction finished in: %.2f %s\n",
  run_time,
  units(run_time)
))

########## Save output and wrap-up ############################################
cat('Saving sample data\n')
# record initial time
start_time <- Sys.time()
# write samples in a csv

cols_loc_time <- sim.obj$quantiles %>%
  filter(time > t1) %>%
  select(lon, lat, time, site_name, actuals.cf, forecast.cf, matches("quant"))

sim.obj$samples %>%
  t() %>%
  as.data.frame() %>%
  setNames(paste0("sim", 1:ncol(.))) %>%
  bind_cols(cols_loc_time) %>%
  # setNames(paste0(
  #   "t_",
  #   seq(t1, t1 + 23 * 3600, by = "hour") %>%
  #     format(., "%Y-%m-%d_%H")
  # )) %>%
  write_parquet(
    .,
    file.path(path.samples, paste0(mod.code, t1, ".parquet"))
  )

end_time <- Sys.time()
run_time <- difftime(end_time, start_time)
cat(sprintf(
  "N = %d samples written in: %.2f %s\n",
  n.samples,
  run_time,
  units(run_time)
))

mend_t <- Sys.time()
run_time <- difftime(mend_t, mstart_t)
cat(sprintf(
  "Whole sample process ended for day: %s in %.2f %s\n",
  format(t1, "%Y-%m-%d"),
  run_time,
  units(run_time)
))
