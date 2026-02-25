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
  geom_sf(data = zoom_bbox, fill = NA, colour = "darkred", linewidth = 1) +
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
    bottom = 0.65,
    right = 1.1,
    top = 0.95
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

scots_wf %>% select(lon, lat) %>% distinct() %>% nrow()

corr_df <- spatial_corr_by_distance_fast(
  scots_wf,
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
sorting_models <- data.frame(
  n = 1:length(ordered.levels),
  code = ordered.levels
)

sorted_list <- sorting_models %>%
  right_join(
    cpo_scores,
    by = c("code" = "root_ofolder")
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

# % hyperparameters

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
data_pit <- read.csv("data_pit_st.csv")

my_palette_0 <- ggsci::pal_lancet()(6)

data_pit %>%
  mutate(model = factor(model, levels = ordered.levels)) %>%
  ggplot() +
  geom_abline(aes(slope = 1, intercept = 0), col = "darkgray") +
  stat_ecdf(aes(pit, col = model), na.rm = TRUE, lwd = 0.8) +
  # facet_wrap(~version)+
  scale_color_manual(values = my_palette_0) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.2, .7),
    legend.background = element_blank(), # Makes background completely transparent
    legend.box.background = element_rect(fill = NA, color = NA) # No border
  ) +
  labs(x = "PIT", y = "ECDF")

ggsave(
  file.path(fig_path, "fig_16.eps"),
  width = 3.5,
  height = 3.5
)
# % Reliability diagrams ####
source("aux_funct_ps.R")
# my_palette
my_palette <- pal_lancet()(6)
# scales::show_col(my_palette)
rel.plot <- plot_reliability(
  global_scores = score.tbl.gb,
  model_list = etar_list,
  my_palette = my_palette,
  show.fig = FALSE,
  code_subset = grepv("AR1", ordered.levels, invert = TRUE)
)
ordered.levels <- c("NPB", "NPN", "NPG", "PN", "PG", "EN", "NEN", "nonpar.")

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
  mutate(label = factor(ofolder, levels = ordered.levels))

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
    legend.position.inside = c(.64, .17),
    plot.title = element_text(size = 10), # Adjust title font size
    axis.text = element_text(size = 8), # Adjust axis text font size
    axis.title = element_text(size = 9), # Adjust axis label font size
    legend.text = element_text(size = 8), # Adjust legend text font size
    legend.title = element_text(size = 8), # Adjust legend title font size
    legend.background = element_blank(), # Makes background completely transparent
    legend.box.background = element_rect(fill = NA, color = NA) # No border
  ) +
  guides(col = guide_legend(ncol = 2)) + # Set legend to have 2 columns
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = my_palette) +
  labs(title = "")
p
ggsave(
  filename = file.path(fig_path, "fig_17.eps"),
  plot = p,
  width = 3.5,
  height = 3.5
)
# % Scores

score.tbl <- readRDS("summaries/spatial_scores.rds")

# score.tbl.hour <- lapply(score.tbl, \(x) x$hour)
score.tbl.day <- lapply(score.tbl, \(x) x$day) %>%
  bind_rows() %>%
  filter(grepl("matern-ar1-etaderiv", model)) %>%
  select(date, model, crps) %>%
  pivot_wider(names_from = model, values_from = crps)

score.tbl.gb <- lapply(score.tbl, \(x) x$global) %>%
  bind_rows()


model_scores <- cpo_scores %>%
  select(
    ofolder,
    mod_prefix
  ) %>%
  left_join(
    score.tbl.gb %>%
      select(model, crps, energy, matches("variogram")) %>%
      mutate(
        energy = energy / 35,
        across(matches("variogram"), ~ . / (24 * 35))
      ),
    by = c("mod_prefix" = "model")
  ) %>%
  select(-mod_prefix)

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
