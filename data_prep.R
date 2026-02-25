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

source("aux_funct_ps.R")
source("aux_funct.R")
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
t1 <- "2024-12-31 23:00:00"
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
  filter(abs(error0) <= 0.2, lat <= 56) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

## sample scottish wind farms
set.seed(0)
sample_wf <- sample(scots_summary$site_name, 35)
sample_wf2 <- scots_summary$site_name
scots_summary %>%
  filter(site_name %in% sample_wf) %>%
  ggplot() +
  geom_sf(data = uk_map, fill = "lightgrey", color = "white") +
  geom_sf(aes(geometry = geometry)) +
  theme_map()
ggsave(
  "fig/scottish_35wfsamp_24_map.pdf",
  width = 3.5,
  height = 5
)

ggplot() +
  geom_sf(data = uk_map, fill = "lightgrey", color = "white") +
  geom_sf(
    aes(geometry = geometry, col = "Train"),
    data = scots_summary %>%
      filter((site_name %in% sample_wf))
  ) +
  geom_sf(
    aes(geometry = geometry, col = "OOS"),
    data = scots_summary %>%
      filter(!(site_name %in% sample_wf))
  ) +
  scale_color_manual(values = pal_aaas()(2) %>% rev()) +
  labs(col = "") +
  theme_map()
ggsave(
  "fig/scottish_wfsamp_24_map.pdf",
  width = 3.5,
  height = 5
)
# plot power curve scatter
scots_wf_filtered <- scots_wf %>%
  filter(abs(error0) <= 0.2) %>%
  filter(site_name %in% sample_wf)

scots_wf_filtered2 <- scots_wf %>%
  filter(abs(error0) <= 0.2) %>%
  filter(site_name %in% sample_wf2)
scots_wf_filtered2 %>%
  ggplot() +
  geom_point(aes(x = ws_h, y = norm_potential), alpha = 0.2)

scots_wf_filtered2 <- scots_wf %>%
  mutate(
    rel_error0 = ifelse(
      norm_power_est0 > 0,
      abs(error0) / norm_power_est0,
      abs(error0)
    )
  ) %>%
  filter(abs(error0) <= 0.25) %>%
  # filter(rel_error0 <= 0.2) %>%
  filter(site_name %in% sample_wf2)
scots_wf_filtered2 %>%
  ggplot() +
  geom_point(aes(x = ws_h, y = norm_potential), alpha = 0.2)

first_time <- min(scots_wf$halfHourEndTime)
scots_wf_filtered <- scots_wf %>%
  # filter(abs(error0) <= 0.2) %>%
  filter(site_name %in% sample_wf) %>%
  rename(time = halfHourEndTime) %>%
  mutate(date = as.Date(time), month = month(date), hour = hour(date)) %>%
  # simple time index
  rename(
    actuals.cf = norm_potential,
    forecast.cf = norm_power_est0,
    ws.w = ws_h
  ) %>%
  mutate(
    fcst_group = inla.group(forecast.cf, n = 40, method = "cut"),
    ws.w_group = inla.group(ws.w, n = 40, method = "cut"),
    err.cf = actuals.cf - forecast.cf,
    fd = c(0, diff(forecast.cf)),
    fd_group = inla.group(fd, n = 40, method = "quantile")
  ) %>%
  group_by(site_name) %>%
  arrange(time) %>%
  mutate(
    t = as.numeric(difftime(time, first_time, units = "hours")),
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
    ),
    power.groupf = cut(
      forecast.cf,
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
  file.path("data", "scottish_wfsamp_24.parquet")
)

arrow::write_parquet(
  scots_wf,
  file.path("data", "scottish_wf_24.parquet")
)

scots_wf_filtered %>%
  pull(time) %>%
  range()

scots_wf_filtered %>%
  pull(site_name) %>%
  unique()


# PC estimation ####
# fit 5 parameter power model per site

A <- 1
B <- -12
C <- 8
D <- 0
G <- 1
x <- seq(0, 30, length.out = 100)
y <- pc5param(x, A, B, C, D, G)
plot(x, y, type = "l")

A <- 1
B <- 12
C <- 22
D <- 0
G <- 1
x <- seq(0, 30, length.out = 100)
y <- pc5param(x, A, B, C, D, G)
plot(x, y, type = "l")

A <- 1
B <- -12
C <- 8
C2 <- 22
D <- 0
G <- 1
x <- seq(0, 30, length.out = 100)
y <- pc6param(x, A, B, C, D, G, C2 = 22, B2 = -B * 3)
plot(x, y, type = "l")

## par transformation ####

source("aux_funct.R")
start <- c(
  A_raw = log(max(y)), # log-scale
  B1_raw = log(12),
  C1 = median(x),
  D = 0,
  G_raw = log(1),
  B2_raw = log(24),
  C2 = quantile(x, 0.999, names = FALSE)
)

ws_vec <- pmax(0, scots_wf_filtered2$ws_h)
power_vec <- scots_wf_filtered2$norm_potential

loss(start, ws_vec, power_vec, p = 1, threshold = 20)

fit <- optim(
  par = start,
  fn = loss,
  x = ws_vec,
  y = power_vec,
  # p = 1.2,
  threshold = 20,
  method = "L-BFGS-B",
  control = list(maxit = 500)
)
fit$convergence
fit$par

solution1 <- raw_to_par(fit$par)

loss(fit$par, ws_vec, power_vec)
loss_cutout(fit$par, ws_vec, power_vec)

get_inflection_points(solution1)
infl_approx(solution1)

data.frame(estimate = solution1) %>%
  write.csv("data/pc_7pars_wfsamp_24b.csv", row.names = TRUE)

x_seq <- seq(0, 35, length.out = 200)
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
plot(x_seq, y_seq)

y_seq2 <- pc7param_r %>%
  do.call(
    c(list(x = x_seq), as.list(fit$par))
  )

scots_wf_filtered2 %>%
  ggplot() +
  geom_point(
    aes(x = ws_h, y = norm_potential, col = "observations"),
    alpha = 0.5
  ) +
  geom_line(
    data = data.frame(x = x_seq, y = y_seq),
    aes(
      x = x,
      y = y,
      col = "fit"
    ),
    # color = "red",
    inherit.aes = FALSE
  ) +
  labs(x = "wind speed", y = "normalised power", col = "") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = blues9[c(7, 4)])
# plot time series of wind speed and power
ggsave(
  "fig/scottish_wfsamp_24_pc-est.png",
  width = 6,
  height = 4
)

p1 <- 1
fit_p <- optim(
  par = start,
  fn = loss,
  x = ws_vec,
  y = power_vec,
  p = p1,
  threshold = 22,
  method = "L-BFGS-B",
  control = list(maxit = 500)
)
fit_p$convergence
fit_p$par

loss(fit_p$par, ws_vec, power_vec, p1, 22)
loss_cutout(fit_p$par, ws_vec, power_vec, p1, 22)

fit_p$par %>%
  raw_to_par() %>%
  get_inflection_points()
fit_p$par %>%
  raw_to_par() %>%
  infl_approx()

y_seq_p <- pc7param_r %>%
  do.call(
    c(list(x = x_seq), as.list(fit_p$par))
  )

scots_wf_filtered2 %>%
  ggplot() +
  geom_point(
    aes(x = ws_h, y = norm_potential, col = "observations"),
    alpha = 0.5
  ) +
  geom_line(
    data = data.frame(x = x_seq, y = y_seq_p),
    aes(
      x = x,
      y = y,
      col = "fit"
    ),
    # color = "red",
    inherit.aes = FALSE
  ) +
  labs(x = "wind speed", y = "normalised power", col = "") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = blues9[c(7, 4)])

# Plot 2 weeks of wind farms time series in sample ####
source("aux_funct_ps.R")

t0 <- as.POSIXct("2024-03-01 00:00:00", tz = "UTC")
h = 24 * 14 # 14 days ahead
# scots_wf_filtered$halfHourEndTime %>% tz()
current_site <- scots_summary$site_name[1]

for (current_site in sample_wf) {
  # print(current_site)
  temp_plot <- plot_nh_ahead(
    data = scots_wf_filtered %>% filter(site_name == current_site),
    mycols = c("actuals.cf", "ws.w"),
    col_labels = c("power", "wind speed"),
    t0 = t0,
    h = 24 * 14,
    show.fig = FALSE,
    xvar = "time",
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

for (current_site in sample_wf) {
  # print(current_site)
  temp_plot <- plot_nh_ahead(
    data = scots_wf_filtered %>% filter(site_name == current_site),
    mycols = c("actuals.cf", "forecast.cf"),
    t0 = t0,
    h = 24 * 14,
    show.fig = FALSE,
    xvar = "time"
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

# model list file ####

model_list_fname <- "model_list_spatial.parquet"

model_list_df <- data.frame(
  id = 1:4
) %>%
  mutate(
    family = "gaussian",
    fderiv = "eta",
    out = "error",
    transformation = "normalised",
    response = "err.cf",
    all_combinations = list(
      c("ws.w_group", "etaderiv"),
      c("ws.w_group", "ar1g", "etaderiv"),
      c("ws.w_group", "matern-ar1"),
      c("ws.w_group", "matern-ar1", "etaderiv")
    )
  ) %>%
  mutate(
    # fderiv = ifelse(id == 3, "fd", fderiv),
    fderiv = ifelse((grepl("etaderiv", all_combinations)), fderiv, "fd")
  )

extra_models <- data.frame(
  id = 5:7
) %>%
  mutate(
    family = "gaussian",
    fderiv = "eta",
    out = "error",
    transformation = "normalised",
    response = "err.cf",
    all_combinations = list(
      c("ws.w_group", "fcst_group", "hour", "matern-ar1", "etaderiv"),
      c("ws.w_group", "fcst_group", "matern-ar1", "etaderiv"),
      c("ws.w_group", "hour", "matern-ar1", "etaderiv")
    )
  ) %>%
  mutate(
    # fderiv = ifelse(id == 3, "fd", fderiv),
    fderiv = ifelse((grepl("etaderiv", all_combinations)), fderiv, "fd")
  )

gauss_models <- bind_rows(model_list_df, extra_models)

beta_models <- bind_rows(model_list_df, extra_models) %>%
  mutate(
    family = "beta",
    out = "observed",
    response = "actuals.cf",
    id = 7 + 1:7
  )


write_parquet(
  model_list_df,
  file.path("data", model_list_fname)
)

write_parquet(
  bind_rows(gauss_models, beta_models),
  file.path("data", "model_list_spatial.parquet")
)
# list.files("../p2_replication/eddie_code/", pattern = "model_list")
# read_rds("../p2_replication/eddie_code/model_listtop.rds")
bind_rows(gauss_models, beta_models)

# list.files("data", "model")
no_winds <- read_parquet("data/model_list_spatial.parquet") %>%
  filter(
    !map_lgl(all_combinations, ~ any(str_detect(.x, "fcst_group")))
  )


no_winds <- no_winds %>%
  mutate(
    all_combinations = map(
      all_combinations,
      ~ gsub("ws.w_group", "fcst_group", .x)
    )
  ) %>%
  # pull(all_combinations) %>%
  # ungroup() %>%
  # View()
  rowwise() %>%
  mutate(
    readable_feat = paste(unlist(all_combinations), collapse = ", ")
  ) %>%
  ungroup() %>%
  mutate(id = 14 + 1:nrow(.))

read_parquet("data/model_list_spatial.parquet") %>%
  mutate(all_combinations = as.list(all_combinations)) %>%
  bind_rows(no_winds) %>%
  rowwise() %>%
  mutate(
    readable_feat = paste(unlist(all_combinations), collapse = ", ")
  ) %>%
  ungroup() %>%
  # View()
  write_parquet(
    file.path("data", "model_list_spatial.parquet")
  )


## New model list v2.0 ####

model_list_fname <- "model_list_spatial_v2.parquet"

model_list_df <- data.frame(
  id = 1:3
) %>%
  mutate(
    family = "gaussian",
    fderiv = "eta",
    out = "error",
    transformation = "normalised",
    response = "err.cf",
    all_combinations = list(
      c("ws.w_group", "fcst_group", "hour", "ar2g", "etaderiv"),
      c("ws.w_group", "fcst_group", "fd_group", "hour", "ar2g"),
      c("ws.w_group", "fcst_group", "hour", "matern-ar1", "etaderiv")
    )
  ) %>%
  mutate(
    # fderiv = ifelse(id == 3, "fd", fderiv),
    fderiv = ifelse((grepl("etaderiv", all_combinations)), fderiv, "fd")
  )

beta_models <- model_list_df %>%
  mutate(
    family = "beta",
    out = "observed",
    response = "actuals.cf",
    id = nrow(model_list_df) + 1:n()
  )

st_beta_variations <- beta_models[rep(3, 3), ] %>%
  mutate(
    all_combinations = list(
      c("ws.w_group", "fcst_group", "matern-ar1", "etaderiv"),
      c("fcst_group", "hour", "matern-ar1", "etaderiv"),
      c("fcst_group", "matern-ar1", "etaderiv")
    ),
    id = 6 + 1:n()
  )

ar_beta_variations <- beta_models[rep(1, 3), ] %>%
  mutate(
    all_combinations = list(
      c("ws.w_group", "fcst_group", "ar2g", "etaderiv"),
      c("fcst_group", "hour", "ar2g", "etaderiv"),
      c("fcst_group", "ar2g", "etaderiv")
    ),
    id = 10 + 1:n()
  )

bind_rows(
  model_list_df,
  beta_models,
  st_beta_variations,
  ar_beta_variations
) %>%
  rowwise() %>%
  mutate(
    readable_feat = paste(unlist(all_combinations), collapse = ", ")
  ) %>%
  ungroup() %>%
  # View()
  write_parquet(
    file.path("data", model_list_fname)
  )
