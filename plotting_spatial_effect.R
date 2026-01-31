source("aux_funct_ps.R")
model_fname <- "r_err.cf_f_gaussian_eta_feat_ws.w_group-matern-ar1-etaderiv.rds"
model_path <- "~/Documents/proj2/spatial/model_objects/"
mod.temp <- readRDS(file.path(model_path, model_fname))

stack_fname <- paste0("misc/stack_", model_fname)
stack <- readRDS(stack_fname)
inla.stack.index(stack, "wf.stack")$data %>% head()
plot.effects.spatial(mod.temp)

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

# full set of indices
idat <- inla.stack.index(stack, 'wf.data')$data
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
  mutate(fitted = mod.temp$summary.fitted.values$mean[idat_pred] + forecast.cf)

fit_obs_df %>%
  ggplot() +
  geom_line(aes(time, fitted, col = "model")) +
  geom_line(aes(time, actuals.cf, col = "observed")) +
  facet_wrap(~site_name) +
  scale_color_manual(
    values = c("darkred", "darkblue"),
    breaks = c("model", "observed")
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
par(mfrow = c(4, 3), mar = c(0, 0, 0, 0))
for (j in 1:24) {
  xmean[[j]][-idx_inside] <- NA
  book.plot.field(
    list(x = projgrid$x, y = projgrid$y, z = xmean[[j]]),
    zlim = round(range(unlist(xmean), na.rm = TRUE), 1)
  )
}
