## Packages ####

require(tidyverse)
require(parallel)
require(ggsci)
require(kableExtra)
require(arrow)
require(sf)
source("aux_funct.R")
source("aux_funct_ps.R")

theme_set(theme_bw())

## data figures/tables ###

input_data <- "data/scottish_wfsamp_24.parquet"
data.scaled <- read_parquet(input_data)

proba <- c(0, 0.2, 0.8, 1)
with(
  data.scaled,
  quantile(actuals.cf, proba, na.rm = TRUE)
)
with(
  data.scaled,
  quantile(forecast.cf, proba, na.rm = TRUE)
)
prob.labels <- c(
  paste0("low", round((proba[2] - proba[1]) * 100), "%"),
  "mid",
  paste0("high", round((proba[4] - proba[3]) * 100), "%")
)

data.scaled <- data.scaled %>%
  mutate(
    power.groupf = cut(
      forecast.cf,
      breaks = quantile(forecast.cf, proba, na.rm = TRUE),
      include.lowest = TRUE,
      labels = prob.labels
    )
  )
prop.table(
  table(data.scaled$power.group, data.scaled$power.groupf),
  margin = 1
)
data_day_summary <- data.scaled %>%
  mutate(date = date(time)) %>%
  group_by(date) %>%
  summarise(
    across(c(actuals.cf, forecast.cf), \(x) mean(x, na.rm = TRUE))
  ) %>%
  mutate(
    power.group = cut(
      actuals.cf,
      breaks = quantile(actuals.cf, proba, na.rm = TRUE),
      include.lowest = TRUE,
      labels = prob.labels
    ),
    power.groupf = cut(
      forecast.cf,
      breaks = quantile(forecast.cf, proba, na.rm = TRUE),
      include.lowest = TRUE,
      labels = prob.labels
    )
  )
tab <- with(
  data_day_summary,
  prop.table(
    table(power.group, power.groupf),
    margin = 1
  )
)

tab_df <- as.data.frame.matrix(tab)

kable(
  tab_df * 100,
  format = "latex",
  digits = 1,
  booktabs = TRUE,
  caption = paste(
    "Row-normalised contingency table comparing power groups defined",
    "from observed capacity factors (power.group) and forecast capacity",
    "factors (power.groupf). Entries show the conditional proportion",
    "of observations assigned to each forecast power group given the",
    "observed power group."
  ),
  col.names = paste("Forecast", colnames(tab_df))
)

## Data section figures ####

### Wind farm map ####
require(rnaturalearth)
require(rnaturalearthdata)
require(ggthemes)

# Wind farm locations
wind_points <- data.scaled %>%
  select(lon, lat) %>%
  distinct() %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# UK polygon
uk <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf") %>%
  filter(admin == "United Kingdom")

myxlim <- c(-6.5, -0.5)
myylim = c(54, 57)
zoom_bbox <- st_as_sfc(
  st_bbox(
    c(
      xmin = myxlim[1],
      xmax = myxlim[2],
      ymin = myylim[1],
      ymax = myylim[2]
    ),
    crs = 4326
  )
)
p_main <- ggplot() +
  geom_sf(data = uk, fill = "grey90", colour = "grey40") +
  geom_sf(data = wind_points, colour = "darkblue") +
  coord_sf(xlim = myxlim, ylim = myylim) +
  theme_map()
p_locator <- ggplot() +
  geom_sf(data = uk, fill = "grey90", colour = "grey40") +
  geom_sf(data = zoom_bbox, fill = NA, colour = "darkred", linewidth = 0.8) +
  theme_void()
library(patchwork)
require(ggspatial)
p_main +
  annotation_scale(
    location = "bl", # bottom-left
    width_hint = 0.3, # fraction of plot width
    unit_category = "metric"
  ) +
  inset_element(
    p_locator,
    left = 0.65,
    bottom = 0.5,
    right = 1,
    top = 1.2
  )
ggsave(
  "fig/fig_5bis.eps",
  width = 3.5,
  height = 3
)


### Correlation ####
# Spatial correlation heatmap
cor_matrix <- data.scaled %>%
  select(time, site_name, actuals.cf) %>%
  pivot_wider(names_from = site_name, values_from = actuals.cf) %>%
  select(-time) %>%
  cor(use = "pairwise.complete.obs")

cor_melted <- reshape2::melt(cor_matrix)

p1 <- ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6)
  ) +
  labs(x = NULL, y = NULL, title = "Spatial correlation of actuals.cf")

# Temporal ACF of cross-sectional mean
mean_cf <- data.scaled %>%
  group_by(time) %>%
  summarise(mean_actuals = mean(err.cf, na.rm = TRUE))

p2 <- forecast::ggAcf(mean_cf$mean_actuals, lag.max = 25) +
  labs(title = "", x = "Lag (hours)") +
  theme_minimal()

p3 <- forecast::ggPacf(mean_cf$mean_actuals, lag.max = 25) +
  labs(title = "", x = "Lag (hours)") +
  theme_minimal()

ggsave(
  "fig/fig_5bis2a.eps",
  p2,
  width = 3.5,
  height = 3
)

ggsave(
  "fig/fig_5bis2b.eps",
  p3,
  width = 3.5,
  height = 3
)

p1 + p2


source("aux_funct.R")

scots_wf <- read_parquet("data/scottish_wf_24.parquet") %>%
  rename(
    actuals.cf = norm_potential,
    err0 = norm_power_est0,
    time = halfHourEndTime
  )

dist_mat <- data.scaled %>%
  select(lon, lat) %>%
  distinct() %>%
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
  select(x, y) %>%
  dist()
dist_mat %>%
  range()
dist_mat %>% as.numeric() %>% density() %>% plot()
dist_mat %>% as.numeric() %>% hist()
scots_wf %>% select(lon, lat) %>% distinct() %>% nrow()

corr_df <- spatial_corr_by_distance_fast(
  data.scaled,
  value_col = "actuals.cf",
  n_bins = 15,
  bin_type = "quantile",
  time_col = "time",
  ribbon_probs = c(0.05, 0.95)
)

p_scot_cor <- ggplot(corr_df, aes(x = dist_mean, y = corr_mean)) +
  geom_ribbon(
    aes(ymin = corr_lower, ymax = corr_upper),
    alpha = 0.5,
    fill = blues9[3]
  ) +
  geom_line(color = blues9[5]) +
  geom_point(color = blues9[7]) +
  # geom_line(aes(y = threshold_95), linetype = "dashed", color = "red") +
  # geom_line(aes(y = -threshold_95), linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Distance (km)",
    y = "Correlation",
    # title = "Spatial correlation vs distance with variability and significance"
  ) +
  theme_minimal()

ggsave("fig/fig_5bis3.pdf", p_scot_cor, width = 3.5, height = 3)


## model figures/table ####

fig_path <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/proj2/prob-scenarios-main-doc/fig_clean/"
main_folder <- "~/Documents/proj2/spatial"
model_path <- "~/Documents/proj2/spatial/model_objects"
list.files(model_path, pattern = "ar2g|matern-ar1")

# mod.temp <- readRDS(
#   file.path(
#     model_path,
#     "r_err.cf_f_gaussian_eta_feat_ws.w_group-ar1g-etaderiv.rds"
#   )
# )
# mod.temp$.args$formula

# model list ####
cpo_files <- c(
  "summaries/top_stmodel_cpo_v2.csv",
  "summaries/top_stmodel_cpo_NPBv2.csv",
  "summaries/top_stmodel_cpo_ARNPBv2.csv"
)


cpo_scores <- read.csv("summaries/top_stmodel_cpo.csv") %>%
  mutate(
    sample_path = file.path(
      main_folder,
      ofolder,
      "sample"
    ),
    mod_prefix = gsub("\\.rds", "_t", mod.file.name)
  ) %>%
  filter(!grepl("AR1", ofolder))

ordered.levels <- c(
  "PR-NPB",
  "AR1-PR-NPB",
  "ST-NPB",
  "ST-PR-NPB",
  "PR-NEN",
  "AR1-PR-NEN",
  "ST-NEN",
  "ST-PR-NEN"
)
cpo_scores <- lapply(
  cpo_files,
  \(f) read.csv(f)
) %>%
  bind_rows() %>%
  mutate(
    sample_path = file.path(
      main_folder,
      ofolder,
      "sample"
    ),
    mod_prefix = gsub("\\.rds", "_t", mod.file.name),
    root_ofolder = sub("_.+$", "", ofolder)
  ) %>%
  group_by(
    root_ofolder
  ) %>%
  filter(
    ifelse(
      !is.na(mean_log_gcpo),
      mean_log_gcpo == min(mean_log_gcpo),
      mean_log_cpo == min(mean_log_cpo)
    )
  )

ordered.levels <- c(
  "AR2-NPB",
  "AR2-PR-NPB",
  "ST-PR-NPB",
  "AR2-NEN",
  "AR2-PR-NEN",
  "ST-PR-NEN"
)
ordered.labels <- sub("AR2-", "", ordered.levels)
sorting_models <- data.frame(
  n = 1:length(ordered.levels),
  code = ordered.levels
)

sorted_list <- sorting_models %>%
  right_join(
    cpo_scores,
    by = c("code" = "root_ofolder")
  )
sorted_list %>%
  write.csv(
    "summaries/sorted_stmodels.csv",
    row.names = FALSE
  )

# % effects ####

fnames <- sorted_list %>%
  mutate(mod_prefix = sub("_t", "", mod_prefix) %>% paste0(., ".rds")) %>%
  pull(mod_prefix)

## ST-PR-NPB ####
# rm(mod.temp)
mod.temp <- readRDS(file.path(model_path, fnames[3]))
get_mod_stats("", "", mod.temp)
effects.list <- names(mod.temp$summary.random)
excluded <- c("t", "eta", "eta.1", "eta.2", "spatial")

t <- lapply(
  effects.list[!effects.list %in% excluded],
  \(effect) plot.effects(mod.temp, effect, show.fig = T)
)

t[[1]]$fig +
  labs(x = "Point forecast", y = "Estimated effect", title = "") +
  scale_x_continuous(
    # n.breaks = 6,
    breaks = pretty(t[[1]]$fig$data$group, n = 4), # use the scaled x range#
    labels = function(x) {
      round(x, 1)
    }
  )
figwidth <- 3
figheight <- 2.5

# ggsave(file.path(fig_path, "fig_13a.eps"), width = figwidth, height = figheight)
ggsave(file.path(fig_path, "fig_13a.pdf"), width = figwidth, height = figheight)


## ST-PR-NEN ####
figheight <- 3
figwidth <- 3.5
mod.temp <- readRDS(file.path(model_path, fnames[6]))
get_mod_stats("", "", mod.temp)
effects.list <- names(mod.temp$summary.random)
excluded <- c("t", "eta", "eta.1", "eta.2", "spatial")

eff_list <- effects.list[!effects.list %in% excluded]
t <- lapply(
  eff_list,
  \(effect) plot.effects(mod.temp, effect, show.fig = T)
)

t[[1]]$fig +
  labs(x = "Wind speed", y = "Estimated effect", title = "") +
  scale_x_continuous(
    # n.breaks = 6,
    breaks = pretty(t[[1]]$fig$data$group, n = 4), # use the scaled x range#
    labels = function(x) {
      round(x, 1)
    } # back-transform
  )

# ggsave(file.path(fig_path, "fig_6d.eps"), width = figwidth, height = figheight)
ggsave(file.path(fig_path, "fig_13c.pdf"), width = figwidth, height = figheight)

t[[2]]$fig +
  labs(x = "Point forecast", y = "Estimated effect", title = "") +
  scale_x_continuous(
    # n.breaks = 6,
    breaks = pretty(t[[2]]$fig$data$group, n = 4), # use the scaled x range#
    labels = function(x) {
      round(x, 1)
    } # back-transform
  )

# ggsave(file.path(fig_path, "fig_6d.eps"), width = figwidth, height = figheight)
ggsave(file.path(fig_path, "fig_13b.pdf"), width = figwidth, height = figheight)


t[[3]]$fig +
  labs(x = "Hour", y = "Estimated effect", title = "") +
  scale_x_continuous(
    # n.breaks = 6,
    breaks = pretty(t[[3]]$fig$data$group, n = 4), # use the scaled x range#
    labels = function(x) {
      round(x, 1)
    } # back-transform
  )

# ggsave(file.path(fig_path, "fig_6d.eps"), width = figwidth, height = figheight)
ggsave(file.path(fig_path, "fig_13d.pdf"), width = figwidth, height = figheight)


## Model Densities ####

# % hyperparameters

model_path <- "~/Documents/proj2/spatial/model_objects"

mod_names <- c(
  "r_actuals.cf_f_beta_eta_feat_fcst_group-matern-ar1-etaderiv.rds",
  "r_err.cf_f_gaussian_eta_feat_ws.w_group-fcst_group-hour-matern-ar1-etaderiv.rds"
)

fname <- "summaries/densities_summaries.rds"
if (!file.exists(fname)) {
  summaries.dens <- lapply(mod_names, \(fname) {
    modeltemp <- readRDS(file.path(model_path, fname))

    c(modeltemp$marginals.fixed, modeltemp$marginals.hyperpar)
  })
  saveRDS(summaries.dens, fname)
} else {
  summaries.dens <- readRDS(fname)
}

source("aux_funct_ps.R")

### NPB ####
densities <- plot.densities_regime(
  sum_dens = summaries.dens[1],
  model_tag = "NPB",
  ncol = 3
)
par_names <- summaries.dens[[1]] %>% names()
par_names_print <- c(
  "Intercept",
  "'Precision Beta'",
  # "Precision~gaussian",
  # "Precision~s[w](.)",
  "Precision~s[f](.)",
  # "Precision~s[m](.)",
  "Range",
  "Spatial~field~sigma",
  "AR~rho",
  "Power~ramp~effect"
)

names(par_names_print) <- par_names
densities +
  labs(title = "", x = "") +
  facet_wrap(
    ~parameter,
    scales = "free",
    ncol = 3,
    labeller = as_labeller(par_names_print, default = label_parsed)
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7), # Adjust axis text font size
    panel.spacing.x = unit(1, "cm") # increase spacing
  ) +
  scale_color_manual(values = "gray25")

ggsave(file.path(fig_path, "fig_S4.eps"), width = 8, height = 4.5)

### NEN ####
source("aux_funct_ps.R")
loglist <- c("Precision for fcst_group", "Precision for ws.w_group")
densities <- plot.densities_regime(
  sum_dens = summaries.dens[2],
  model_tag = "NEN",
  ncol = 3,
  loglist = loglist
)
par_names <- summaries.dens[[2]] %>% names()
par_names_print <- c(
  "Intercept",
  # "'Precision Beta'",
  "Precision~Gaussian",
  "Precision~s[w](.)",
  "Precision~s[f](.)",
  "Precision~s[h](.)",
  "Range",
  "Spatial~field~sigma",
  "AR~rho",
  "Power~ramp~effect"
)
names(par_names_print) <- par_names

# densities %>% str(levels = 1)
densities +
  labs(title = "", x = "") +
  facet_wrap(
    ~parameter,
    scales = "free",
    ncol = 3,
    labeller = as_labeller(par_names_print, default = label_parsed)
  ) +
  facetted_pos_scales(
    x = list(
      parameter %in% loglist ~ scale_x_log10()
    )
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7), # Adjust axis text font size
    panel.spacing.x = unit(1, "cm") # increase spacing
  ) +
  scale_color_manual(values = "gray25")

ggsave(file.path(fig_path, "fig_S5.eps"), width = 8, height = 4.5)


hyper_names <- mod.temp$summary.hyperpar %>% rownames()
hyper_names_print <- c(
  "'Precision Beta'",
  "Precision~s[f](.)",
  "Range",
  "Spatial~field~sigma",
  "rho",
  "Power~ramp~effect"
)
# sapply(hyper_names_print, latex2exp::TeX)
names(hyper_names_print) <- hyper_names
phyper <- plot.hyper.dens(mod.temp, show.fig = FALSE)
phyper +
  labs(title = "", x = "") +
  facet_wrap(
    ~parameter,
    scales = "free",
    ncol = 2,
    labeller = as_labeller(hyper_names_print, default = label_parsed)
  )

ggsave(
  file.path(fig_path, "fig_14a.pdf"),
  width = figwidth * 2,
  height = figheight * 2
)

# % hyperparameters

hyper_names <- mod.temp$summary.hyperpar %>% rownames()
hyper_names_print <- c(
  "Precision~gaussian",
  "Precision~s[w](.)",
  "Precision~s[f](.)",
  "Precision~s[h](.)",
  "Range",
  "Spatial~field~sigma",
  "rho",
  "Power~ramp~effect"
)
# sapply(hyper_names_print, latex2exp::TeX)
names(hyper_names_print) <- hyper_names
phyper <- plot.hyper.dens(mod.temp, show.fig = FALSE)
phyper +
  labs(title = "", x = "") +
  facet_wrap(
    ~parameter,
    scales = "free",
    ncol = 2,
    labeller = as_labeller(hyper_names_print, default = label_parsed)
  )

ggsave(
  file.path(fig_path, "fig_14b.pdf"),
  width = figwidth * 2,
  height = figheight * 2
)

# % scenarios 1st day ####

ofolder <- sorted_list$ofolder
date <- "2024-07-05"
plot_type <- "sim.small"
mod_prefix <- gsub("\\.rds", "_t", sorted_list$mod.file.name)

sampfigs <- mapply(
  \(folder, prefix, i) {
    fname <- file.path(
      "/home/s2441782/Documents/proj2/spatial",
      folder,
      "fig",
      paste0(prefix, date, " 23:00:00", "_", plot_type, ".png")
    )
    # print(fname)
    # file.exists(fname) %>% print()
    # fname
    knitr::include_graphics(fname)
    figname <- file.path(
      fig_path,
      paste0("fig_15", i, ".png")
    )
    if (file.exists(fname)) {
      file.copy(from = fname, to = figname, overwrite = TRUE)
    } else {
      warning(paste("File not found:", fname))
    }
    # ggsave(
    #   figname,
    #   width = figwidth,
    #   height = figheight
    # )
    invisible(figname)
  },
  folder = ofolder,
  prefix = mod_prefix,
  i = c("a", "b", "c", "d", "e", "f")
)
# sampfigs[1] %>% unname()
# "/home/s2441782/Documents/proj2/spatial/ST-PR-NPB/fig/r_actuals.cf_f_beta_eta_feat_fcst_group-matern-ar1-etaderiv_t2024-07-05_sim.small.png"
# fname <- file.path("~/Documents/proj2/spatial", )

# % scenarios 2nd day

# % PIT diagrams ####

fnames <- sorted_list %>%
  mutate(mod_prefix = sub("_t", "", mod_prefix) %>% paste0(., ".rds")) %>%
  pull(mod_prefix)

data_pit <- lapply(
  fnames,
  \(file) {
    # browser()
    mod.temp <- readRDS(file.path(model_path, file))

    preffect_present <- any(grepl("^eta", mod.temp$summary.random %>% names()))

    if (preffect_present) {
      n_eta <- length(mod.temp$cpo$pit) / 2
      pit_vals <- mod.temp$cpo$pit[-c(1:n_eta)]
    } else {
      n_eta <- length(mod.temp$cpo$pit)
      pit_vals <- mod.temp$cpo$pit
    }
    rm(mod.temp)
    data.frame(
      id = 1:n_eta,
      pit = pit_vals,
      model = sorted_list$code[which(fnames == file)]
    )
  }
) %>%
  bind_rows()

write.csv(data_pit, "summaries/data_pit_st.csv", row.names = FALSE)
data_pit <- read.csv("summaries/data_pit_st.csv")

my_palette_0 <- ggsci::pal_lancet()(4)

p <- data_pit %>%
  mutate(model = factor(model, levels = ordered.levels)) %>%
  filter(!grepl("AR2-NPB|AR2-NEN", model)) %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0), col = "darkgray") +
  stat_ecdf(aes(pit, col = model), na.rm = TRUE, lwd = 0.8) +
  # facet_wrap(~version)+
  scale_color_manual(
    values = my_palette_0,
    labels = c("PR-NPB", "ST-PR-NPB", "PR-NEN", "ST-PR-NEN")
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.75, .25),
    plot.title = element_text(size = 10), # Adjust title font size
    axis.text = element_text(size = 8), # Adjust axis text font size
    axis.title = element_text(size = 9), # Adjust axis label font size
    legend.text = element_text(size = 8), # Adjust legend text font size
    legend.title = element_text(size = 8), # Adjust legend title font size
    legend.background = element_blank(), # Makes background completely transparent
    legend.box.background = element_rect(fill = NA, color = NA) # No border
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "PIT", y = "ECDF", col = "")
p
ggsave(
  file.path(fig_path, "fig_16.eps"),
  p,
  width = 3,
  height = 3
)
# % Reliability diagrams ####
source("aux_funct_ps.R")
# my_palette
my_palette <- pal_lancet()(6)
# scales::show_col(my_palette)
# rel.plot <- plot_reliability(
#   global_scores = score.tbl.gb,
#   model_list = etar_list,
#   my_palette = my_palette,
#   show.fig = FALSE,
#   code_subset = grepv("AR1", ordered.levels, invert = TRUE)
# )
# ordered.levels <- c("NPB", "NPN", "NPG", "PN", "PG", "EN", "NEN", "nonpar.")

plot.data <- score.tbl.gb %>%
  select(1:29) %>%
  # names()
  pivot_longer(
    cols = coverage05:coverage95,
    names_to = "level",
    values_to = "empirical"
  ) %>%
  mutate(nominal = (substr(level, 9, 10) %>% as.numeric()) / 100) %>%
  right_join(
    cpo_scores %>% select(mod_prefix, ofolder),
    by = c("model" = "mod_prefix")
  ) %>%
  filter(!grepl("AR2-NPB|AR2-NEN", root_ofolder)) %>%
  mutate(
    label = factor(
      root_ofolder,
      levels = ordered.levels,
      labels = ordered.labels
    )
  )

p <- plot.data %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0), col = "darkgray") +
  geom_line(
    aes(
      nominal,
      empirical,
      col = label
    ),
    lwd = 0.8
  ) +
  labs(
    col = "",
    title = "Reliability diagrams",
    x = "Nominal coverage",
    y = "Empirical coverage"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.75, .25),
    plot.title = element_text(size = 10), # Adjust title font size
    axis.text = element_text(size = 8), # Adjust axis text font size
    axis.title = element_text(size = 9), # Adjust axis label font size
    legend.text = element_text(size = 8), # Adjust legend text font size
    legend.title = element_text(size = 8), # Adjust legend title font size
    legend.background = element_blank(), # Makes background completely transparent
    legend.box.background = element_rect(fill = NA, color = NA) # No border
  ) +
  guides(col = guide_legend(ncol = 1)) + # Set legend to have 2 columns
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = my_palette) +
  labs(title = "")
p
ggsave(
  filename = file.path(fig_path, "fig_17.eps"),
  plot = p,
  width = 3,
  height = 3
)
# % Scores

score.tbl <- readRDS("summaries/spatial_scores.rds")


# score.tbl.hour <- lapply(score.tbl, \(x) x$hour)
score.tbl.day <- lapply(score.tbl, \(x) x$day) %>%
  bind_rows() %>%
  filter(grepl("etaderiv", model)) %>%
  # filter(grepl("matern-ar1-etaderiv", model)) %>%
  select(date, model, crps) %>%
  pivot_wider(names_from = model, values_from = crps)

score.tbl.site <- lapply(score.tbl, \(x) x$site) %>%
  bind_rows() %>%
  filter(grepl("etaderiv", model)) %>%
  # filter(grepl("matern-ar1-etaderiv", model)) %>%
  select(site_name, model, crps) %>%
  pivot_wider(names_from = model, values_from = crps)

score.tbl.gb <- lapply(score.tbl, \(x) x$global) %>%
  bind_rows()


cpo_scores <- read.csv("summaries/top_stmodel_cpo.csv") %>%
  mutate(
    sample_path = file.path(
      main_folder,
      ofolder,
      "sample"
    ),
    mod_prefix = gsub("\\.rds", "_t", mod.file.name)
  ) %>%
  filter(!grepl("AR1", ofolder))

ordered.levels <- c(
  "PR-NPB",
  # "AR1-PR-NPB",
  # "ST-NPB",
  "ST-PR-NPB",
  "PR-NEN",
  # "AR1-PR-NEN",
  # "ST-NEN",
  "ST-PR-NEN"
)

# cpo_scores$root_ofolder
model_scores <- cpo_scores %>%
  ungroup() %>%
  select(
    ofolder,
    mod_prefix
  ) %>%
  filter(grepl("etaderiv", mod_prefix)) %>%
  left_join(
    score.tbl.gb %>%
      select(model, crps, energy, matches("variogram")) %>%
      mutate(
        energy = energy / 35,
        across(matches("variogram"), ~ . / (24 * 35))
      ),
    by = c("mod_prefix" = "model")
  ) %>%
  select(-c(mod_prefix))

sorting_models <- data.frame(
  n = 1:length(ordered.levels),
  code = ordered.levels
)

scores1.tbl <- sorting_models %>%
  left_join(model_scores, by = c("code" = "ofolder")) %>%
  select(-n) %>%
  # Apply bolding to all numeric columns
  mutate(across(
    where(is.numeric),
    ~ {
      # 1. Round first so the comparison and display are consistent
      val <- round(.x, 3)
      # 2. Apply bolding to the rounded minimum
      cell_spec(val, "latex", bold = val == min(val, na.rm = TRUE))
    }
  ))

scores1.tbl %>%
  kable(digits = 3, format = "latex", escape = FALSE) %>%
  kable_styling(full_width = F) %>%
  cat()

score.tbl <- readRDS("summaries/spatial_scores_v2.rds")
score.tbl.gb <- lapply(score.tbl, \(x) x$global) %>%
  bind_rows()

# Scenarios Groupped ####

sorted_list <- read.csv(
  "summaries/sorted_stmodels.csv"
)

ofolder <- sorted_list$ofolder
date0 <- "2024-08-05"

plot_type <- "sim.small"
mod_prefix <- gsub("\\.rds", "_t", sorted_list$mod.file.name)

sample_path <- "~/Documents/proj2/spatial"
extension <- " 23:00:00.parquet"

h0 <- 23
t0 <- as.Date(date0) + hours(h0)
fcst_times <- seq.POSIXt(
  from = as.Date(date0) + hours(h0 + 1),
  to = as.Date(date0) + hours(h0 + 24),
  by = "hour"
)

train_data <- history_window(
  data.scaled,
  as.Date(date0) + hours(h0),
  # as.Date("2024-06-30") + hours(h0),
  window = 7,
  units = "days",
  mask = FALSE
) %>%
  arrange(site_name, time)
fcst_ind <- which(train_data$time %in% fcst_times)

cols_replacement <- train_data %>%
  select(time, site_name, actuals.cf, forecast.cf) %>%
  slice(fcst_ind) %>%
  mutate(
    site_name = gsub(
      "wind farm|Wind farm|Wind Farm|Windfarm|\\(Brockloch Rig Phase 2\\)|, Wester Dod Community |Clyde Wind Farm Extension|\\(|\\)| and Dalquhandy Renewable Energy Project| Extension",
      "",
      site_name
    ) %>%
      trimws()
  )

n.samp.large <- 20
set.seed(0)
selected_sim <- sample(1:1000, n.samp.large)
single_wf <- "Fallago Rig"
set_3wf <- c("Aikengall II", "Blackcraig", "Sanquhar Community")
myextension <- "eps"
k <- 6

lapply(
  c(2, 3, 5, 6),
  \(k) {
    samp_temp <- read_parquet(
      file.path(
        sorted_list$sample_path,
        sprintf("%s%s%s", mod_prefix, date0, extension)
      )[k]
    )

    if (k %in% c(1:6)) {
      samp_temp <- samp_temp %>%
        select(-c(site_name, time, actuals.cf, forecast.cf)) %>%
        bind_cols(
          cols_replacement
        )
    }

    pattern_sel <- "^sim"
    # pattern_sel <- "^quant"

    aggr_df <- samp_temp %>%
      select(
        site_name,
        time,
        actuals.cf,
        forecast.cf,
        all_of(selected_sim),
        matches("^quant")
      ) %>%
      group_by(time) %>%
      summarise(
        across(c(actuals.cf, forecast.cf, matches(pattern_sel)), \(x) {
          mean(x, na.rm = TRUE)
        })
      ) %>%
      mutate(
        across(
          matches(pattern_sel),
          \(x) if (grepl("err.cf", mod_prefix[k])) x + forecast.cf else x
        )
      ) %>%
      pivot_longer(
        cols = c(actuals.cf, forecast.cf, matches(pattern_sel)),
        names_to = "type",
        values_to = "value"
      ) %>%
      mutate(
        value = pmin(1, pmax(0, value)),
        type = factor(
          type,
          levels = c(
            paste0("sim", selected_sim),
            "forecast.cf",
            "actuals.cf"
          )
        )
      )
    aggr_df %>%
      ggplot() +
      geom_line(aes(time, value, col = type)) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_color_manual(
        values = (c("darkred", "darkblue", rep("grey50", 30))),
        breaks = c("actuals.cf", "forecast.cf", paste0("sim", selected_sim[1])),
        labels = c("observed", "forecast", "scenarios")
      ) +
      labs(col = "", y = "") +
      theme(
        legend.position = ifelse(k == 6, "inside", "none"),
        legend.position.inside = c(0.25, 0.75),
        legend.background = element_blank(), # Makes background completely transparent
        legend.box.background = element_rect(fill = NA, color = NA) # No border
      ) +
      scale_x_datetime(date_labels = "%H:%M")

    ggsave(
      sprintf("fig_fcst/aggr_%s%s.%s", sorted_list$code[k], date0, myextension),
      width = 3.5,
      height = 3
    )

    selectwf_df <- samp_temp %>%
      select(
        site_name,
        time,
        actuals.cf,
        forecast.cf,
        all_of(selected_sim),
        matches("^quant")
      ) %>%
      filter(site_name %in% set_3wf) %>%
      mutate(
        across(
          matches(pattern_sel),
          \(x) if (grepl("err.cf", mod_prefix[k])) x + forecast.cf else x
        )
      ) %>%
      pivot_longer(
        cols = c(actuals.cf, forecast.cf, matches(pattern_sel)),
        names_to = "type",
        values_to = "value"
      ) %>%
      mutate(
        value = pmin(1, pmax(0, value)),
        type = factor(
          type,
          levels = c(
            paste0("sim", selected_sim),
            "forecast.cf",
            "actuals.cf"
          )
        )
      )
    bind_rows(
      aggr_df %>% mutate(site_name = "Aggregation"),
      selectwf_df
    ) %>%
      ggplot() +
      geom_line(aes(time, value, col = type)) +
      coord_cartesian(ylim = c(0, 1)) +
      facet_wrap(~site_name, ncol = 4) +
      scale_color_manual(
        values = (c("darkred", "darkblue", rep("grey50", 30))),
        breaks = c("actuals.cf", "forecast.cf", paste0("sim", selected_sim[1])),
        labels = c("observed", "forecast", "scenarios")
      ) +
      labs(col = "", y = "Wind Generation % of Capacity") +
      theme(
        legend.position = ifelse(k == 6, "inside", "none"),
        legend.position.inside = c(0.07, 0.8),
        legend.background = element_blank(), # Makes background completely transparent
        legend.box.background = element_rect(fill = NA, color = NA) # No border
      ) +
      scale_x_datetime(date_labels = "%H")

    ggsave(
      sprintf(
        "fig_fcst/selectwf_waggr_%s%s.%s",
        sorted_list$code[k],
        date0,
        myextension
      ),
      width = 10,
      height = 3
    )

    samp_temp %>%
      select(
        site_name,
        time,
        actuals.cf,
        forecast.cf,
        all_of(selected_sim),
        matches("^quant")
      ) %>%
      mutate(
        across(
          matches(pattern_sel),
          \(x) if (grepl("err.cf", mod_prefix[k])) x + forecast.cf else x
        )
      ) %>%
      pivot_longer(
        cols = c(actuals.cf, forecast.cf, matches(pattern_sel)),
        names_to = "type",
        values_to = "value"
      ) %>%
      mutate(
        value = pmin(1, pmax(0, value)),
        type = factor(
          type,
          levels = c(
            paste0("sim", selected_sim),
            "forecast.cf",
            "actuals.cf"
          )
        )
      ) %>%
      ggplot() +
      geom_line(aes(time, value, col = type)) +
      coord_cartesian(ylim = c(0, 1)) +
      facet_wrap(~site_name) +
      scale_color_manual(
        values = (c("darkred", "darkblue", rep("grey50", 30))),
        breaks = c("actuals.cf", "forecast.cf", paste0("sim", selected_sim[1])),
        labels = c("observed", "forecast", "scenarios")
      ) +
      labs(col = "", y = "Wind Generation % of Capacity") +
      theme(
        legend.position = ifelse(k == 6, "inside", "none"),
        legend.position.inside = c(0.9, 0.1),
        legend.background = element_blank(), # Makes background completely transparent
        legend.box.background = element_rect(fill = NA, color = NA) # No border
      ) +
      scale_x_datetime(date_labels = "%H")

    ggsave(
      sprintf("fig_fcst/35wf_%s%s.%s", sorted_list$code[k], date0, myextension),
      width = 10,
      height = 9
    )
  }
)

# samp_temp <- read_parquet(
#   file.path(
#     sorted_list$sample_path,
#     sprintf("%s%s%s", mod_prefix, date0, extension)
#   )[k]
# )

# samp_temp <- samp_temp %>%
#   select(-c(site_name, time, actuals.cf, forecast.cf)) %>%
#   bind_cols(
#     cols_replacement
#   )

# n.samp.large <- 20
# set.seed(1)
# selected_sim <- sample(1:1000, n.samp.large)
# # samp_temp %>%
# #   select(site_name, time, actuals.cf, all_of(selected_sim)) %>%
# #   head()

# pattern_sel <- "^sim"
# # pattern_sel <- "^quant"

# samp_temp %>%
#   select(
#     site_name,
#     time,
#     actuals.cf,
#     forecast.cf,
#     all_of(selected_sim),
#     matches("^quant")
#   ) %>%
#   group_by(time) %>%
#   summarise(
#     across(c(actuals.cf, forecast.cf, matches(pattern_sel)), \(x) {
#       mean(x, na.rm = TRUE)
#     })
#   ) %>%
#   mutate(
#     across(
#       matches(pattern_sel),
#       \(x) if (grepl("err.cf", mod_prefix[k])) x + forecast.cf else x
#     )
#   ) %>%
#   pivot_longer(
#     cols = c(actuals.cf, forecast.cf, matches(pattern_sel)),
#     names_to = "type",
#     values_to = "value"
#   ) %>%
#   mutate(value = pmin(1, pmax(0, value))) %>%
#   ggplot() +
#   geom_line(aes(time, value, col = type)) +
#   coord_cartesian(ylim = c(0, 1)) +
#   # facet_wrap(~site_name) +
#   scale_color_manual(
#     values = c("actuals.cf" = "darkred", "forecast.cf" = "darkblue")
#   ) +
#   labs(col = "") +
#   theme(legend.position = "bottom")

# debug(plot_actuals_model)
# test <- plot_actuals_model(
#   data = data.scaled,
#   samp_temp[, 1:1000] %>% t(),
#   response = "actuals.cf",
#   t1 = t0,
#   h = 24,
#   resp.lab = "Wind Generation",
#   plot.type = "sim",
#   n.sim.plot = n.samp.large,
#   show.fig = TRUE,
#   legend.opt = "forecast",
#   clipping = TRUE,
#   spatial = TRUE,
#   fcst_points = fcst_ind
# )

# samp_temp

# plotting spatial latent field using vignette ####
require(INLA)
require(fmesher)
require(ggthemes)
require(inlabru)

# Wind farm locations
wind_points <- data.scaled %>%
  select(lon, lat) %>%
  distinct() %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# UK polygon
uk <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf") %>%
  filter(admin == "United Kingdom")
sorted_list <- read.csv(
  "summaries/sorted_stmodels.csv"
)

ofolder <- sorted_list$ofolder
date0 <- "2024-08-05"
t1 <- "2024-08-05 23:00:00" %>% as.POSIXct(tz = "UTC")
mod_prefix <- gsub("\\.rds", "_t", sorted_list$mod.file.name)

fig_path <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/proj2/prob-scenarios-main-doc/fig_clean/"
main_folder <- "~/Documents/proj2/spatial"
model_path <- "~/Documents/proj2/spatial/model_objects"

mod.temp <- readRDS(file.path(model_path, sorted_list$mod.file.name[6]))

## ----samples-------------------------------------------------------------
nn <- 20
s <- inla.posterior.sample(
  n = nn,
  mod.temp,
  intern = TRUE,
  seed = 0,
  add.names = FALSE
)

## Find the values of latent field "i" in samples from mesh1
contents <- mod.temp$misc$configs$contents
effect <- "spatial"
id.effect <- which(contents$tag == effect)
ind.effect <- contents$start[id.effect] - 1 + (1:contents$length[id.effect])


# Obtain predictions at the nodes of mesh2
spde <- mod.temp$.args$data$wf.spde
mesh <- mod.temp$.args$data$wf.spde$mesh
loc1 = mesh$loc[, 1:2]
loc2 = mesh$loc[, 1:2]
n = mesh$n

mtch = match(data.frame(t(loc2)), data.frame(t(loc1)))
idx.c = which(!is.na(mtch))
idx.u = setdiff(1:mesh$n, idx.c)
p = c(idx.u, idx.c)

ypred.mesh2 = matrix(c(NA), mesh$n, nn)

m <- n - length(idx.c)
iperm <- numeric(m)


## ----meanrf--------------------------------------------------------------
# cor(fit_obs_df$actuals.cf, fit_obs_df$fitted)

## ----projgrid------------------------------------------------------------
stepsize <- 2

data_masked <- history_window(
  data.scaled,
  t1,
  window = 7,
  units = "days",
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
  st_drop_geometry()
# coords <- data_masked %>% select(x, y) %>% unique() %>% as.matrix()

coords <- data.scaled %>%
  distinct(lon, lat) %>%
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
  select(x, y) %>%
  as.matrix()
# bnd <- fm_extensions(coords, convex = c(-.10, -.15))

loc_unique <- data_masked %>%
  distinct(x, y) %>%
  as.matrix()
bnd <- fm_extensions(loc_unique, convex = c(-.2, -.20))
bnd[[1]] %>% st_coordinates()
# bnd <- fm_extensions(coords, convex = c(-.1, -.20))
ggplot() + geom_sf(data = bnd[[1]])
wf.spde <- mod.temp$.args$data$wf.spde
mesh <- wf.spde$mesh
# bnd <- fmesher::fm_segm(mesh, boundary = TRUE)
# require(fmesher)
# bnd <- inla.mesh.boundary(mesh)
# ggplot() + gg(data = mesh)
nxy <- round(c(diff(range(coords[, 1])), diff(range(coords[, 2]))) / stepsize)
# projgrid <- inla.mesh.projector(
#   mesh,
#   xlim = range(coords[, 1]),
#   ylim = range(coords[, 2]),
#   dims = nxy
# )
projgrid <- fm_evaluator(
  mesh,
  xlim = range(coords[, 1]),
  ylim = range(coords[, 2]),
  dims = nxy
)

projgrid <- fm_evaluator(
  mesh,
  xlim = range(st_coordinates(bnd[[1]])[, 1]),
  ylim = range(st_coordinates(bnd[[1]])[, 2]),
  dims = nxy
)
# expand.grid(projgrid$x, projgrid$y) %>%
#   plot()
## ----projpmean-----------------------------------------------------------
time_grid <- seq(
  from = min(data_masked$time),
  to = max(data_masked$time),
  by = "1 hour"
)

data_masked$time_idx <- match(data_masked$time, time_grid)
k <- data_masked$time_idx %>% unique() %>% length()
st.group <- mod.temp$.args$data$st.group[!is.na(mod.temp$.args$data$st.group)]
summary(st.group)

spde_idx <- inla.spde.make.index(
  name = "spatial",
  n.spde = wf.spde$n.spde,
  n.group = length(unique(st.group))
)

# (193 - 23:0) %>% length()
# xfit <- list()
# for (j in (1:24)) {
#   xfit[[j]] <- fm_evaluate(
#     projgrid,
#     mod.temp$summary.fitted.values$mean[spde_idx$spatial.group == k - 24 + j]
#   )
# }
# plot(xfit[1])
# mod.temp$summary.fitted.values$mean %>% length()
# mod.temp$summary.random$spatial$mean %>% length()

xmean <- list()
for (j in (1:24)) {
  xmean[[j]] <- fm_evaluate(
    projgrid,
    mod.temp$summary.random$spatial$mean[spde_idx$spatial.group == k - 24 + j]
  )
}
xmode <- list()
for (j in (1:24)) {
  xmode[[j]] <- fm_evaluate(
    projgrid,
    mod.temp$summary.random$spatial$mode[spde_idx$spatial.group == k - 24 + j]
  )
}

xsd <- list()
for (j in (1:24)) {
  xsd[[j]] <- fm_evaluate(
    projgrid,
    mod.temp$summary.random$spatial$sd[spde_idx$spatial.group == k - 24 + j]
  )
}

xlow <- list()
for (j in (1:24)) {
  xlow[[j]] <- fm_evaluate(
    projgrid,
    mod.temp$summary.random$spatial$`0.025quant`[
      spde_idx$spatial.group == k - 24 + j
    ]
  )
}

xhigh <- list()
for (j in (1:24)) {
  xhigh[[j]] <- fm_evaluate(
    projgrid,
    mod.temp$summary.random$spatial$`0.975quant`[
      spde_idx$spatial.group == k - 24 + j
    ]
  )
}

## ----inout---------------------------------------------------------------
# library(splancs)
# projgrid$lattice$loc %>% plot()
# xy.in <- inout(projgrid$lattice$loc, cbind(PRborder[, 1], PRborder[, 2]))
pts_sf <- st_as_sf(
  data.frame(
    x = projgrid$lattice$loc[, 1],
    y = projgrid$lattice$loc[, 2]
  ),
  coords = c("x", "y"),
  crs = st_crs(bnd)
)

# coords <- bnd[[1]]
# coords_closed <- rbind(coords, coords[1, ])
# bnd_poly <- st_sfc(
#   st_polygon(list(coords_closed)),
#   crs = st_crs(pts_sf) # ensure same CRS
# )
# plot(coords)
# plot(coords, pch = 16)
# text(coords, labels = seq_len(nrow(coords)), pos = 3, cex = 0.7)
# plot(coords_closed)

# plot(mesh)
# plot(bnd[[2]])
# plot(pts_sf)
# plot(bnd[[1]])
# mesh$loc %>% str()
# plot(bnd[[1]], col = "blue", lwd = 2)
# lines(mesh$loc[bnd[[1]][[2]], ], col = "red", lwd = 2)
inside <- st_within(pts_sf, bnd[[1]], sparse = FALSE)[, 1]
# inside <- st_within(pts_sf, bnd_poly, sparse = FALSE)[, 1]
idx_inside <- which(inside)
times <- seq(1, 24, 2)

pdf(
  "fig/speff_ST-PR-NEN.pdf",
  width = 10,
  height = 7
)
par(mfrow = c(4, 3), mar = c(1, 1, 1, 2))

# times <- c(1, 12, 24)
for (j in times) {
  xmean[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xmean[[j]]),
    zlim = round(range(unlist(xmean), na.rm = TRUE), 1),
    main = sprintf(
      "Time: %s",
      format(as.POSIXct("2024-07-01", tz = "UTC") + hours(j), "%H:%M")
    )
  )
  # plot(mesh, add = TRUE)
  points(loc_unique)
}
dev.off()


pdf(
  "fig/speffmode_ST-PR-NEN.pdf",
  width = 10,
  height = 7
)
par(mfrow = c(4, 3), mar = c(1, 1, 1, 2))


for (j in times) {
  xmode[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xmode[[j]]),
    zlim = round(range(unlist(xmode), na.rm = TRUE), 1),
    main = sprintf(
      "Time: %s",
      format(as.POSIXct("2024-07-01", tz = "UTC") + hours(j), "%H:%M")
    )
  )
  # plot(mesh, add = TRUE)
  points(loc_unique)
}
dev.off()

pdf(
  "fig/speffsd_ST-PR-NEN.pdf",
  width = 10,
  height = 7
)
par(mfrow = c(4, 3), mar = c(1, 1, 1, 2))

for (j in times) {
  xsd[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xsd[[j]]),
    zlim = round(range(unlist(xsd), na.rm = TRUE), 1),
    main = sprintf(
      "Time: %s",
      format(as.POSIXct("2024-07-01", tz = "UTC") + hours(j), "%H:%M")
    ),
    # mesh = mesh,
    col = book.color.c2()
  )
  # plot(mesh, add = TRUE)
  points(loc_unique)
}
dev.off()

pdf(
  "fig/spefflow_ST-PR-NEN.pdf",
  width = 10,
  height = 7
)
par(mfrow = c(4, 3), mar = c(1, 1, 1, 2))

for (j in times) {
  xlow[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xlow[[j]]),
    zlim = round(range(unlist(xlow), na.rm = TRUE), 1),
    main = sprintf(
      "Time: %s",
      format(as.POSIXct("2024-07-01", tz = "UTC") + hours(j), "%H:%M")
    ),
    # mesh = mesh,
    col = book.color.c()
  )
  # plot(mesh, add = TRUE)
  points(loc_unique)
}
dev.off()

pdf(
  "fig/speffhigh_ST-PR-NEN.pdf",
  width = 10,
  height = 7
)
par(mfrow = c(4, 3), mar = c(1, 1, 1, 2))

for (j in times) {
  xhigh[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xhigh[[j]]),
    zlim = round(range(unlist(xhigh), na.rm = TRUE), 1),
    main = sprintf(
      "Time: %s",
      format(as.POSIXct("2024-07-01", tz = "UTC") + hours(j), "%H:%M")
    ),
    # mesh = mesh,
    col = book.color.c()
  )
  # plot(mesh, add = TRUE)
  points(loc_unique)
}
dev.off()


## Summary plot with mean field, sd, and quantiles ####
times <- c(1, 12, 24)

col_lims <- round(range(c(unlist(xlow), unlist(xhigh)), na.rm = TRUE), 1)

pdf(
  "fig/sp_summary_ST-PR-NEN.pdf",
  width = 10,
  height = 7
)
par(mfrow = c(4, 3), mar = c(1, 1, 1, 2))
for (j in times) {
  xmean[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xmean[[j]]),
    zlim = col_lims,
    main = sprintf(
      "Mean , h=%s",
      format(as.POSIXct("2024-07-01", tz = "UTC") + hours(j), "%H:%M")
    )
  )
  points(loc_unique)
}
for (j in times) {
  xsd[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xsd[[j]]),
    zlim = round(range(unlist(xsd), na.rm = TRUE), 1),
    main = sprintf(
      "Std. dev. %s",
      format(as.POSIXct("2024-07-01", tz = "UTC") + hours(j), "%H:%M")
    ),
    col = book.color.c2()
  )
  points(loc_unique)
}
for (j in times) {
  xlow[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xlow[[j]]),
    zlim = col_lims,
    main = sprintf(
      "Lower quant. %s",
      format(as.POSIXct("2024-07-01", tz = "UTC") + hours(j), "%H:%M")
    )
  )
  points(loc_unique)
}
for (j in times) {
  xhigh[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xhigh[[j]]),
    zlim = col_lims,
    main = sprintf(
      "Upper quant. %s",
      format(as.POSIXct("2024-07-01", tz = "UTC") + hours(j), "%H:%M")
    ),
    # mesh = mesh,
    col = book.color.c()
  )
  points(loc_unique)
}
dev.off()

grid <- expand.grid(
  x = projgrid$x,
  y = projgrid$y
)

# 2) flatten mask to match the same ordering
mask_vec <- as.vector(idx_inside) # must be 111 x 86

df_list <- list()

for (j in times) {
  # flatten fields
  v_mean <- as.vector(xmean[[j]])
  v_sd <- as.vector(xsd[[j]])
  v_low <- as.vector(xlow[[j]])
  v_high <- as.vector(xhigh[[j]])

  # apply mask
  v_mean[!mask_vec] <- NA
  v_sd[!mask_vec] <- NA
  v_low[!mask_vec] <- NA
  v_high[!mask_vec] <- NA

  df_list[[paste0("mean_", j)]] <- data.frame(
    x = grid$x,
    y = grid$y,
    value = v_mean,
    stat = "Posterior mean",
    lead = paste0("+", j, " h")
  )

  df_list[[paste0("sd_", j)]] <- data.frame(
    x = grid$x,
    y = grid$y,
    value = v_sd,
    stat = "Posterior std. dev.",
    lead = paste0("+", j, " h")
  )

  df_list[[paste0("low_", j)]] <- data.frame(
    x = grid$x,
    y = grid$y,
    value = v_low,
    stat = "2.5% quantile",
    lead = paste0("+", j, " h")
  )

  df_list[[paste0("high_", j)]] <- data.frame(
    x = grid$x,
    y = grid$y,
    value = v_high,
    stat = "97.5% quantile",
    lead = paste0("+", j, " h")
  )
}

df <- bind_rows(df_list)

# enforce ordering for facets
df$stat <- factor(
  df$stat,
  levels = c(
    "Posterior mean",
    "Posterior std. dev.",
    "2.5% quantile",
    "97.5% quantile"
  )
)

myxlim <- c(-6, -0)
myylim <- c(54.5, 56.5)

uk_crop <- uk %>%
  st_crop(
    xmin = myxlim[1],
    xmax = myxlim[2],
    ymin = myylim[1],
    ymax = myylim[2]
  )
uk_27700 <- st_transform(uk_crop, 27700)

# ggplot(uk_crop) +
#   geom_sf() +
#   coord_sf(xlim = myxlim, ylim = myylim)
df <- df %>%
  mutate(
    x = x * 1000,
    y = y * 1000
  )
loc_df <- as.data.frame(loc_unique) %>%
  setNames(c("x", "y")) %>%
  mutate(
    x = x * 1000,
    y = y * 1000
  )
df$lead <- factor(df$lead, levels = paste0("+", times, " h"))
p <- ggplot(df, aes(x = x, y = y, fill = value)) +
  # geom_sf(
  #   data = uk_27700,
  #   fill = "grey90",
  #   colour = "grey40",
  #   inherit.aes = FALSE
  # ) +
  geom_raster(
    # alpha = 0.8
  ) +
  facet_grid(stat ~ lead, switch = "y") +
  # coord_equal() +
  # scale_fill_viridis_c(limits = col_lims, na.value = NA) +
  scale_fill_viridis_c(limits = col_lims, na.value = NA, option = "A") +
  theme_map() +
  theme(
    strip.text = element_text(size = 10),
    axis.title = element_blank(),
    legend.position = "right",
    legend.position.inside = c(-0.1, 0.1)
  ) +
  labs(fill = "")
p
p <- p +
  geom_point(
    data = as.data.frame(loc_df),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    size = 0.4
  )
ggsave("fig/sp_summary_ST-PR-NEN.pdf", p, width = 10, height = 7)


p <- ggplot() +
  geom_sf(data = uk_27700, fill = "grey90", colour = "grey40") +
  geom_raster(data = df, aes(x = x, y = y, fill = value)) +
  geom_point(
    data = loc_df,
    aes(x = x, y = y),
    size = 0.4
  ) +
  facet_grid(stat ~ lead, switch = "y") +
  coord_sf(crs = 27700, xlim = myxlim, ylim = myylim) +

  scale_fill_viridis_c(limits = col_lims, na.value = NA) +
  theme_map() +
  theme(
    strip.text = element_text(size = 10),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  labs(fill = "")


uk_crop_3857 <- st_transform(uk_crop, 3857)
df_sf <- st_as_sf(df, coords = c("x", "y"), crs = 27700) %>%
  st_transform(3857) %>%
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  )

p <- ggplot() +

  # geom_sf(
  #   data = uk_crop_3857,
  #   fill = "grey90",
  #   colour = "grey40",
  #   inherit.aes = FALSE
  # ) +

  geom_tile(
    data = df_sf,
    aes(x = x, y = y, fill = value)
  ) +

  facet_grid(stat ~ lead, switch = "y") +

  scale_fill_viridis_c(limits = col_lims) +

  theme_map() +
  theme(
    strip.text = element_text(size = 10),
    axis.title = element_blank(),
    legend.position = "right"
  ) +

  labs(fill = "")
p
