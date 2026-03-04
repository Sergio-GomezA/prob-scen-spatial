score.tbl.day %>% View()
scen.tbl %>% pull(actuals.cf)[1]
abline(v = (scen.tbl %>% pull(actuals.cf))[1])


checks %>%
  group_by(site_name) %>%
  summarise(
    across(c(crps, energy, matches("variogram")), \(x) mean(x, na.rm = T))
  ) %>%
  view()

checks$time %>% unique()
checks %>%
  filter(site_name %in% c("Whiteside Hill")) %>%
  mutate(hour = hour(time)) %>%
  group_by(hour) %>%
  summarise(
    across(c(crps, energy, matches("variogram")), \(x) mean(x, na.rm = T))
  ) %>%
  view()
