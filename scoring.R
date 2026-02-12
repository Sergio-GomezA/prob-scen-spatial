# list of models
require(data.table)
require(tidyverse)
require(parallel)

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

main_folder <- "~/Documents/proj2/spatial"
cpo_scores <- read.csv("summaries/top_stmodel_cpo.csv") %>%
  mutate(
    sample_path = file.path(
      main_folder,
      ofolder,
      "sample"
    ),
    mod_prefix = gsub("\\.rds", "_t", mod.file.name)
  )

source("aux_funct_ps.R")
# debug(model.scoring.reliability)
# debug(stats_hour_model)
# debug(functions_to_scenarios)

# data
input_data <- "data/scottish_wfsamp_24.parquet"
data.scaled <- read_parquet(input_data) %>%
  mutate(forecast.cf_orig = forecast.cf, forecast = forecast.cf * capacity)

# first day to forecast
t0 <- as.POSIXct("2024-07-01", tz = "UTC")

# days to run
time_seq <- seq(
  from = t0,
  to = t0 + 60 * 24 * 3600,
  by = "day"
)
require(arrow)
test <- model.scoring.reliability(
  model.name = cpo_scores$mod_prefix[1],
  time_seq = time_seq,
  sample.path = cpo_scores$sample_path[1],
  compressed = TRUE,
  ext = ".parquet"
)

# hour
score.tbl <- mcmapply(
  \(model, path) {
    tryCatch(
      model.scoring.reliability(
        model,
        time_seq,
        path,
        compressed = TRUE,
        ext = ".parquet"
      ),
      error = function(e) {
        message("Error when scoring samples: ", e$message)
        return(NULL)
      }
    )
  },
  model = cpo_scores$mod_prefix,
  path = cpo_scores$sample_path,
  mc.cores = available_cores() - 2,
  SIMPLIFY = FALSE
)

saveRDS(score.tbl, "summaries/spatial_scores.rds")
