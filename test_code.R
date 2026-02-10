# pcov <- myspde.posterior(mod.temp, "spatial", "matern.covariance")
# # pcov %>% names()
# pcov %>%
#   ggplot() +
#   geom_ribbon(
#     aes(x, ymin = `q0.025%`, ymax = `q0.975%`),
#     alpha = 0.5
#   ) +
#   gg(pcov)

# pcor <- myspde.posterior(mod.temp, "spatial", "matern.correlation")
# pcor %>%
#   ggplot() +
#   geom_ribbon(
#     aes(x, ymin = `q0.025%`, ymax = `q0.975%`),
#     alpha = 0.5
#   ) +
#   gg(pcor)

model_fnames <- list.files("~/Documents/proj2/spatial/model_objects")

test_df <- data.frame(
  mod.file.name = model_fnames
) %>%
  mutate(
    ofolder = case_when(
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
  )

scores_path <- "~/Documents/proj2/spatial/etaderiv"
cpo_scores <- list.files(
  scores_path
) %>%
  lapply(., \(x) fread(file.path(scores_path, x)) %>% mutate(fname = x)) %>%
  bind_rows() %>%
  mutate(
    mod.file.name = paste0(gsub("var_scores_|\\.csv", "", fname), ".rds")
  ) %>%
  right_join(., test_df, by = "mod.file.name") %>%
  select(ofolder, mod.file.name, everything()) %>%
  group_by(
    ofolder
  ) %>%
  filter(
    !grepl("ws.w_group-fcst_group", mod.file.name)
  ) %>%
  filter(
    waic == min(waic)
  )

write.csv(cpo_scores, "summaries/top_stmodel_cpo.csv", row.names = F)

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
  TRUE ~ "Other"
)
