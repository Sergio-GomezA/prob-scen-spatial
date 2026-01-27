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

source("aux_funct.R")
source("aux_funct_ps.R")
source("functions_probscen.R")

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
ref_catalog_2025 <- fread(
  file.path("data/ref_catalog_wind_2025_era.csv.gz")
)
source("aux_funct.R")
# inla.setOption(num.threads = paste0(n.cores, ":1"))

# plot location of wind farms
pwr_curv_df %>% head()

t0 <- "2024-01-01 00:00:00"
t1 <- "2024-06-30 23:00:00"
scots_wf <- pwr_curv_df %>%
  left_join(
    ref_catalog_2025 %>%
      select(bmUnit, country) %>%
      unique(),
    by = "bmUnit"
  ) %>%
  filter(
    country == "Scotland",
    halfHourEndTime >= t0,
    halfHourEndTime <= t1
  ) %>%
  group_by(lon, lat, halfHourEndTime) %>%
  summarise(
    site_name = first(site_name),
    tech_typ = first(tech_typ),
    # ws_h = mean(ws_h),
    across(c(ws_h, ws10, ws100, wd10, wd100), mean),
    across(c(potential, power_est0, capacity), sum),
    .groups = "drop"
  ) %>%
  mutate(
    norm_potential = pmin(1, potential / capacity),
    norm_power_est0 = power_est0 / capacity,
    error0 = norm_potential - norm_power_est0
  ) %>%
  mutate(
    site_id = as.integer(factor(site_name)),
    pos_val = ifelse(
      norm_potential > 0 & norm_potential < 1,
      pmin(pmax(norm_potential, 1e-6), 1 - 1e-6),
      NA_real_
    ),
    is_zero = as.integer(norm_potential == 0),
    generic_logit = qlogis(
      pmin(pmax(norm_power_est0, 1e-6), 1 - 1e-6)
    )
  )
scots_wf$site_name %>% unique() %>% length()

uk_map <- ne_countries(
  scale = "medium",
  country = "United Kingdom",
  returnclass = "sf"
)
# plot locations
scots_summary <- scots_wf %>%
  filter(tech_typ == "Wind Onshore") %>%
  group_by(site_name, lon, lat) %>%
  summarise(
    n = n(),
    across(c(ws_h, matches("ws|wd|potential|power_est0|norm_|err")), mean)
  ) %>%
  filter(abs(error0) <= 0.2) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

scots_summary %>%
  ggplot() +
  geom_sf(data = uk_map, fill = "lightgrey", color = "white") +
  geom_sf(aes(geometry = geometry))

# plot power curve scatter
scots_wf_filtered <- scots_wf %>%
  filter(abs(error0) <= 0.2) %>%
  filter(site_name %in% scots_summary$site_name)
scots_wf_filtered %>%
  ggplot() +
  geom_point(aes(x = ws_h, y = norm_potential), alpha = 0.1)


scots_wf_filtered <- scots_wf %>%
  # filter(abs(error0) <= 0.2) %>%
  filter(site_name %in% scots_summary$site_name) %>%
  rename(time = halfHourEndTime) %>%
  mutate(date = as.Date(time)) %>%
  # simple time index
  rename(actuals.cf = norm_potential, forecast.cf = power_est0, ws.w = ws_h) %>%
  group_by(site_name) %>%
  arrange(time) %>%
  mutate(
    t = as.numeric(difftime(time, first(time), units = "hours")),
    fcst_group = inla.group(forecast.cf, n = 20, method = "cut"),
    ws.w_group = inla.group(ws.w, n = 20, method = "cut"),
    err.cf = actuals.cf - forecast.cf,
    fd = c(0, diff(forecast.cf)),
    fd_group = inla.group(fd, n = 10, method = "quantile")
  ) %>%
  ungroup()

# grouping days for adverse situations
proba <- c(0, 0.2, 0.8, 1)
with(
  scots_wf_filtered,
  quantile(actuals.cf, proba, na.rm = TRUE)
)

prob.labels <- c(
  paste0("low", round((proba[2] - proba[1]) * 100), "%"),
  "mid",
  paste0("high", round((proba[4] - proba[3]) * 100), "%")
)

grouping.days <- scots_wf_filtered %>%
  # filter(time >= t1) %>%
  group_by(date) %>%
  summarise(
    err.cf = mean(err.cf, na.rm = TRUE),
    # day.range = diff(range(actuals.cf, na.rm = TRUE)),
    actuals.cf = mean(actuals.cf, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Dynamically set labels
    err.group = cut(
      err.cf,
      breaks = quantile(err.cf, proba, na.rm = TRUE),
      include.lowest = TRUE,
      labels = prob.labels
    ),
    # day.range.group = cut(
    #   day.range,
    #   breaks = quantile(day.range, proba, na.rm = TRUE),
    #   include.lowest = TRUE,
    #   labels = prob.labels
    # ),
    power.group = cut(
      actuals.cf,
      breaks = quantile(actuals.cf, proba, na.rm = TRUE),
      include.lowest = TRUE,
      labels = prob.labels
    )
  )

scots_wf_filtered <- scots_wf_filtered %>%
  left_join(
    grouping.days %>% rename(mean.err = err.cf, mean.pow = actuals.cf),
    by = "date"
  )

arrow::write_parquet(
  scots_wf_filtered,
  file.path("data", "scotish_wf_24.parquet")
)

scots_wf_filtered %>%
  pull(halfHourEndTime) %>%
  range()

scots_wf_filtered %>%
  pull(site_name) %>%
  unique()


# fit 5 parameter power model per site

pc5param <- function(x, A, B, C, D, G) {
  D + (A - D) / (1 + (x / C)^B)^G
}

A <- 1
B <- -12
C <- 8
D <- 0
G <- 1
x <- seq(0, 30, length.out = 100)
y <- pc5param(x, A, B, C, D, G)
plot(x, y, type = "l")

pc5param <- function(x, A, B, C, D, G) {
  D + (A - D) / (1 + (x / C)^B)^G
}

A <- 1
B <- 12
C <- 22
D <- 0
G <- 1
x <- seq(0, 30, length.out = 100)
y <- pc5param(x, A, B, C, D, G)
plot(x, y, type = "l")

pc6param <- function(x, A, B, C, D, G, C2, B2 = -B) {
  ifelse(
    x < (C + C2) / 2,
    D + (A - D) / (1 + (x / C)^B)^G,
    D + (A - D) / (1 + (x / C2)^B2)^G
  )
}

A <- 1
B <- -12
C <- 8
C2 <- 22
D <- 0
G <- 1
x <- seq(0, 30, length.out = 100)
y <- pc6param(x, A, B, C, D, G, C2 = 22, B2 = -B * 3)
plot(x, y, type = "l")


pc8param <- function(x, A, B1, C1, D, G1, C2, B2, G2) {
  D + (A - D) / (1 + (x / C1)^B1)^G / (1 + (x / C2)^B2)^G2
}
A <- 1
B <- -12
B2 <- 24
C <- 8
C2 <- 22
D <- 0
G <- 1
G2 <- 1
x <- seq(0, 30, length.out = 100)
y <- pc8param(x, A, B, C, D = 0, G, C2, B2, G2)
plot(x, y, type = "l")


pc8param <- function(x, A, B1, C1, D, G1, C2, B2, G2) {
  D + (A - D) / (1 + (x / C1)^B1)^G / (1 + (x / C2)^B2)^G2
}
loss <- function(par, x, y) {
  # browser()
  y_hat <- pc8param(
    x,
    A = par["A"],
    B = par["B"],
    C = par["C"],
    D = par["D"],
    G = par["G"],
    C2 = par["C2"],
    B2 = par["B2"],
    G2 = par["G2"]
  )
  sum((y - y_hat)^2)
}

x <- pmax(1e-6, scots_wf_filtered$ws_h)
y <- scots_wf_filtered$norm_potential
start <- c(
  A = max(y),
  C = median(x),
  B = -12,
  D = 0,
  G = 1,
  C2 = quantile(x, 0.999, names = FALSE),
  B2 = 24,
  G2 = 1
)

fit <- optim(
  par = start,
  fn = loss,
  x = x,
  y = y,
  method = "L-BFGS-B"
)

#### par transformation
pc7param_r <- function(x, A_raw, B1_raw, C1, D, G_raw, C2, B2_raw) {
  B1 <- -exp(B1_raw) # always negative
  B2 <- exp(B2_raw) # always positive
  A <- D + exp(A_raw) # ensures A > D
  G <- exp(G_raw) # ensures G > 0

  D +
    (A - D) /
      (1 + (x / C1)^B1)^G /
      (1 + (x / C2)^B2)^G
}
raw_to_par <- function(par) {
  B1 <- -exp(par["B1_raw"]) # always negative
  B2 <- exp(par["B2_raw"]) # always positive
  A <- par["D"] + exp(par["A_raw"]) # ensures A > D
  G <- exp(par["G_raw"]) # ensures G > 0
  c(
    A = A,
    B1 = B1,
    C1 = par["C1"],
    D = par["D"],
    G = G,
    C2 = par["C2"],
    B2 = B2
  )
}
loss <- function(par, x, y) {
  # browser()
  y_hat <- pc7param_r(
    x,
    A = par["A_raw"],
    B1_raw = par["B1_raw"],
    C1 = par["C1"],
    D = par["D"],
    G_raw = par["G_raw"],
    C2 = par["C2"],
    B2_raw = par["B2_raw"]
  )
  sum((y - y_hat)^2)
}
start <- c(
  A_raw = log(max(y)), # log-scale
  B1_raw = log(12),
  C1 = median(x),
  D = 0,
  G_raw = log(1),
  C2 = quantile(x, 0.999, names = FALSE),
  B2_raw = log(24)
)
loss(start, x, y)

fit <- optim(
  par = start,
  fn = loss,
  x = x,
  y = y,
  method = "L-BFGS-B"
)
fit$par
raw_to_par <- function(par) {
  B1 <- -exp(unname(par["B1_raw"]))
  B2 <- exp(unname(par["B2_raw"]))
  A <- unname(par["D"]) + exp(unname(par["A_raw"]))
  G <- exp(unname(par["G_raw"]))

  c(
    A = A,
    B1 = B1,
    C1 = unname(par["C1"]),
    D = unname(par["D"]),
    G = G,
    C2 = unname(par["C2"]),
    B2 = B2
  )
}
raw_to_par(fit$par)

x_seq <- seq(0, 30, length.out = 200)
y_seq <- pc7param_r(
  x_seq,
  A_raw = fit$par["A_raw"],
  B1_raw = fit$par["B1_raw"],
  C1 = fit$par["C1"],
  D = fit$par["D"],
  G_raw = fit$par["G_raw"],
  C2 = fit$par["C2"],
  B2_raw = fit$par["B2_raw"]
)

scots_wf_filtered %>%
  ggplot() +
  geom_point(aes(x = ws_h, y = norm_potential), alpha = 0.1) +
  geom_line(
    data = data.frame(x = x_seq, y = y_seq),
    aes(
      x = x,
      y = y
    ),
    color = "red",
    inherit.aes = FALSE
  )
# plot time series of wind speed and power

source("aux_funct_ps.R")

t0 <- as.POSIXct("2024-03-01 00:00:00", tz = "UTC")
h = 24 * 14 # 14 days ahead
# scots_wf_filtered$halfHourEndTime %>% tz()
current_site <- scots_summary$site_name[1]

for (current_site in scots_summary$site_name) {
  # print(current_site)
  temp_plot <- plot_nh_ahead(
    data = scots_wf_filtered %>% filter(site_name == current_site),
    mycols = c("norm_potential", "ws_h"),
    col_labels = c("power", "wind speed"),
    t0 = t0,
    h = 24 * 14,
    show.fig = FALSE,
    xvar = "halfHourEndTime",
    axis_2 = TRUE
  )

  ggsave(
    filename = sprintf(
      "fig_ts/ts_wf-%s_t0-%s_h-%03d.png",
      current_site,
      format(t0, "%Y%m%d%H%M"),
      h
    ),
    plot = temp_plot$plot +
      # scale_colour_discrete(labels = c("Observed", "Wind speed")) +
      coord_cartesian(ylim = c(0, NA)) +
      labs(x = "time", y = "Normalised power") +
      scale_x_datetime(),
    width = 8,
    height = 6
  )
}

# get aux functions from last time

# plot time series for a few wind farms
t0 <- as.POSIXct("2024-03-01 00:00:00", tz = "UTC")
h = 24 * 14 # 14 days ahead
# scots_wf_filtered$halfHourEndTime %>% tz()
current_site <- scots_summary$site_name[1]

for (current_site in scots_summary$site_name) {
  # print(current_site)
  temp_plot <- plot_nh_ahead(
    data = scots_wf_filtered %>% filter(site_name == current_site),
    mycols = c("norm_potential", "norm_power_est0"),
    t0 = t0,
    h = 24 * 14,
    show.fig = FALSE,
    xvar = "halfHourEndTime"
  )

  ggsave(
    filename = sprintf(
      "fig_ts/ts_wf-%s_t0-%s_h-%03d.png",
      current_site,
      format(t0, "%Y%m%d%H%M"),
      h
    ),
    plot = temp_plot$plot +
      labs(x = "time", y = "Normalised power") +
      scale_x_datetime(),
    width = 8,
    height = 6
  )
}
# plot PC derived estimate and observed power

# recover wind speed for times without power observations
