# read model

# effects

# hyperparameters

# samples

data_masked %>%
  # group_by(site_id)
  ggplot(aes(x = time, y = fd, group = site_id, col = site_id)) +
  geom_line() +
  facet_wrap(~site_name)


data_test <- data_masked %>%
  mutate(site_id = as.factor(site_id)) %>%
  group_by(site_id) %>%
  mutate(
    p_diff1 = lag(actuals.cf) - lag(actuals.cf, n = 2, default = NA),
    p_diff3 = lag(actuals.cf) - lag(actuals.cf, n = 4, default = NA)
  ) %>%
  ungroup()

data_test %>%
  ggplot(aes(
    x = p_diff1,
    y = actuals.cf,
    group = site_id,
    col = power.group
  )) +
  geom_point(alpha = 0.3) +
  facet_wrap(~site_name)

data_test %>%
  ggplot(aes(
    x = p_diff3,
    y = actuals.cf,
    group = site_id,
    col = power.group
  )) +
  geom_point(alpha = 0.3) +
  facet_wrap(~site_name)


mod1 <- lm(actuals.cf ~ p_diff1 * site_id, data = data_test)
summary(mod1)


mod2 <- lm(actuals.cf ~ p_diff3 * site_id, data = data_test)
summary(mod2)
