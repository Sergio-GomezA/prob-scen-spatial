pcov <- myspde.posterior(mod.temp, "spatial", "matern.covariance")
# pcov %>% names()
pcov %>%
  ggplot() +
  geom_ribbon(
    aes(x, ymin = `q0.025%`, ymax = `q0.975%`),
    alpha = 0.5
  ) +
  gg(pcov)

pcor <- myspde.posterior(mod.temp, "spatial", "matern.correlation")
pcor %>%
  ggplot() +
  geom_ribbon(
    aes(x, ymin = `q0.025%`, ymax = `q0.975%`),
    alpha = 0.5
  ) +
  gg(pcor)
