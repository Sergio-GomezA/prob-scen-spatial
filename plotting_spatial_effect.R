require(tidyverse)
require(INLA)
theme_set(theme_bw())

source("aux_funct_ps.R")
model_fname <- "r_err.cf_f_gaussian_eta_feat_ws.w_group-matern-ar1-etaderiv.rds"
model_path <- "~/Documents/proj2/spatial/model_objects/"
mod.temp <- readRDS(file.path(model_path, model_fname))

stack_fname <- paste0("misc/stack_", model_fname)
stack <- readRDS(stack_fname)
# inla.stack.index(stack, "wf.stack")$data %>% head()
# plot.effects.spatial(mod.temp)

stack %>% attributes()
stack$A %>% dim()

data_masked %>% nrow() * 2
mod.temp$.args$data$wf.spde$n.spde
mod.temp$.args$data$wf.spde$mesh$loc %>% head()

data_masked$time_idx %>% unique() %>% length()
mod.temp$.args$data$eta %>% length()
mod.temp$.args$data$st.group %>% unique() %>% length()
# inla.spde.make.A()

(7 * 24 + 25) * mod.temp$.args$data$wf.spde$n.spde
224 * 193 / 2
mod.temp$summary.fitted.values %>% row.names() %>% length()
grepv("APredictor", row.names(mod.temp$summary.fitted.values)) %>% length()

mod.temp$.args$data$Y_err.cf %>% length()
names_fitted <- row.names(mod.temp$summary.fitted.values)
names_fitted[!grepl("APredictor", names_fitted)] %>% length()

substr(names_fitted, 8, 12) %>% unique()
substr(row.names(mod.temp$summary.fitted.values), 1, 5) %>% unique()


## ----rfidx---------------------------------------------------------------
t1 <- "2024-06-30 23:00:00" %>% as.POSIXct(tz = "UTC")
input_data <- "data/scottish_wfsamp_24.parquet"
data.scaled <- read_parquet(input_data)

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
# full set of indices
# idat <- inla.stack.index(stack, 'wf.data')$data
idat <- stack$data$index$wf.data
n <- length(idat) / 2
# excluding fakezeros
idat_resp <- idat[n + 1:n]

# next day indices
resp_Y <- mod.temp$.args$data$Y_err.cf[n + 1:n, 2]
idat_pred <- which(is.na(resp_Y)) + n
fit_obs_df <- pred_df <- data_masked %>%
  slice(idat_pred - n) %>%
  select(time, site_name) %>%
  left_join(data.scaled, by = c("time", "site_name")) %>%
  mutate(
    fitted = mod.temp$summary.fitted.values$mean[idat_pred] + forecast.cf,
    linpred = mod.temp$summary.linear.predictor$mean[idat_pred] + forecast.cf
  )

fit_obs_df %>%
  ggplot() +
  geom_line(aes(time, fitted, col = "model")) +
  # geom_line(aes(time, linpred, col = "lin pred")) +
  geom_line(aes(time, forecast.cf, col = "point forecast")) +
  geom_line(aes(time, actuals.cf, col = "observed")) +
  facet_wrap(~site_name) +
  scale_color_manual(
    values = c("darkred", "darkblue", "gray50"),
    breaks = c("observed", "model", "point forecast")
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.92, 0.05),
    legend.background = element_blank(), # Makes background completely transparent
    legend.box.background = element_rect(fill = NA, color = NA), # No border
    axis.text = element_text(angle = 90)
  ) +
  labs(col = "", x = "hour", y = "normalised power") +
  scale_x_datetime(date_labels = "%H:%M")
ggsave(
  "fig/ts_wfsamp_24-07-01.pdf",
  width = 10,
  height = 7
)
mod.temp$summary.fitted.values$mean %>% length()
mod.temp$summary.fitted.values[idat_pred, ] %>% head()
mod.temp$summary.fitted.values %>% row.names() %>% substr(., 1, 10) %>% unique()

grepv("fitted\\.APred", row.names(mod.temp$summary.fitted.values)) %>%
  length() /
  2
require(INLA)
test_samples <- windpow.samples <- inla.posterior.sample(
  n = 100,
  result = mod.temp,
  selection = list(Predictor = idat_pred),
  # num.threads = "1:1",
  seed = 0
)
linpredictor.samples <- sapply(test_samples, \(x) x$latent) %>% t()
precision.samples2 <- sapply(test_samples, \(x) x$hyperpar) %>% t()

precision.samples <- inla.hyperpar.sample(n = 100, result = mod.temp)
phi.samples <- precision.samples2[, 1]
set.seed(0)
post.pred.samples <- sapply(
  1:100,
  \(x) {
    rnorm(
      # add gaussian noise
      ncol(linpredictor.samples),
      linpredictor.samples[x, ],
      1 / sqrt(phi.samples[x])
    )
  }
) %>%
  t()

quants <- post.pred.samples %>%
  apply(., 2, quantile, probs = c(0.025, 0.975))
fit_obs_df %>%
  mutate(
    samp1 = post.pred.samples[1, ] + forecast.cf,
    q_l = quants[1, ] + forecast.cf,
    q_h = quants[2, ] + forecast.cf
  ) %>%
  ggplot() +
  geom_ribbon(
    aes(time, ymin = q_l, ymax = q_h, fill = "95% cred")
  ) +
  geom_line(aes(time, fitted, col = "model")) +
  # geom_line(aes(time, linpred, col = "lin pred")) +
  # geom_line(aes(time, samp1, col = "samp1")) +
  geom_line(aes(time, forecast.cf, col = "point forecast")) +
  geom_line(aes(time, actuals.cf, col = "observed")) +
  facet_wrap(~site_name) +
  scale_color_manual(
    values = c("darkred", "darkblue", "gray50", "gray25"),
    breaks = c("observed", "model", "point forecast", "samp1")
  ) +
  scale_fill_manual(,
    values = blues9[3]
  ) +
  theme(
    legend.position = "bottom",
    legend.position.inside = c(0.92, 0.05),
    legend.background = element_blank(), # Makes background completely transparent
    legend.box.background = element_rect(fill = NA, color = NA), # No border
    axis.text = element_text(angle = 90)
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(col = "", x = "hour", y = "normalised power") +
  scale_x_datetime(date_labels = "%H:%M")
# test_samples[[1]]$latent %>% head()
ggsave(
  "fig/ts2_wfsamp_24-07-01.pdf",
  width = 10,
  height = 10
)

## ----meanrf--------------------------------------------------------------
cor(fit_obs_df$actuals.cf, fit_obs_df$fitted)

## ----projgrid------------------------------------------------------------
stepsize <- 2
# coords <- data_masked %>% select(x, y) %>% unique() %>% as.matrix()

coords <- data_masked %>%
  distinct(x, y) %>%
  as.matrix()
bnd <- fm_extensions(coords, convex = c(-.10, -.15))

data <- mod.temp$.args$data
wf.spde <- data$wf.spde
mesh <- wf.spde$mesh
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
# expand.grid(projgrid$x, projgrid$y) %>%
#   plot()
## ----projpmean-----------------------------------------------------------
k <- data_masked$time_idx %>% unique() %>% length()
st.group <- data$st.group[!is.na(data$st.group)]
summary(st.group)

spde_idx <- inla.spde.make.index(
  name = "spatial",
  n.spde = wf.spde$n.spde,
  n.group = length(unique(st.group))
)

k - 24
# (193 - 23:0) %>% length()
xmean <- list()
for (j in (1:24)) {
  xmean[[j]] <- inla.mesh.project(
    projgrid,
    mod.temp$summary.random$spatial$mean[spde_idx$spatial.group == k - 24 + j]
  )
}

## ----inout---------------------------------------------------------------
library(splancs)
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
inside <- st_within(pts_sf, bnd[[1]], sparse = FALSE)[, 1]
idx_inside <- which(inside)
pdf(
  "fig/speff_ws.w_group_mod_r_err.cf_f_gaussian_eta_feat_ws.w_group-matern-ar1-etaderiv.pdf",
  width = 10,
  height = 7
)
par(mfrow = c(4, 3), mar = c(1, 1, 1, 2))


for (j in seq(2, 24, 2)) {
  xmean[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xmean[[j]]),
    zlim = round(range(unlist(xmean), na.rm = TRUE), 1),
    main = sprintf(
      "Time: %s",
      format(as.POSIXct("2024-07-01", tz = "UTC") + hours(j), "%H:%M")
    )
  )
}
dev.off()
