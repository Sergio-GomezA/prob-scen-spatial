## Packages ####
require(tidyverse)
require(parallel)
require(ggsci)
require(kableExtra)

fig_path <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/proj2/prob-scenarios-main-doc/fig_clean/"
model_path <- "~/Documents/proj2/spatial/model_objects"
# model list ####
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
sorting_models <- data.frame(
  n = 1:length(ordered.levels),
  code = ordered.levels
)

# % effects ####

fnames <- cpo_scores %>%
  mutate(mod_prefix = sub("_t", "", mod_prefix) %>% paste0(., ".rds")) %>%
  pull(mod_prefix)

## ST-PR-NPB ####
mod.temp <- readRDS(file.path(model_path, fnames[1]))
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
mod.temp <- readRDS(file.path(model_path, fnames[4]))
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
    } # back-transform
  )

# ggsave(file.path(fig_path, "fig_6d.eps"), width = figwidth, height = figheight)
ggsave(file.path(fig_path, "fig_13b.pdf"), width = figwidth, height = figheight)

# % hyperparameters

hyper_names <- mod.temp$summary.hyperpar %>% rownames()
hyper_names_print <- c(
  "Precision~gaussian",
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
  file.path(fig_path, "fig_14b.pdf"),
  width = figwidth * 2,
  height = figheight * 2
)

# % scenarios 1st day ####

sorted_list <- sorting_models %>%
  right_join(
    cpo_scores,
    by = c("code" = "ofolder")
  )

ofolder <- sorted_list$code
date <- "2024-07-05"
plot_type <- "sim.small"
mod_prefix <- gsub("\\.rds", "_t", sorted_list$mod.file.name)

sampfigs <- mapply(
  \(folder, prefix, i) {
    fname <- file.path(
      "/home/s2441782/Documents/proj2/spatial",
      folder,
      "fig",
      paste0(prefix, date, "_", plot_type, ".png")
    )
    # print(fname)
    # file.exists(fname)
    # fname
    # knitr::include_graphics(fname)
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
sampfigs[1] %>% unname()
"/home/s2441782/Documents/proj2/spatial/ST-PR-NPB/fig/r_actuals.cf_f_beta_eta_feat_fcst_group-matern-ar1-etaderiv_t2024-07-05_sim.small.png"
fname <- file.path("~/Documents/proj2/spatial", )


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

    data.frame(
      id = 1:n_eta,
      pit = pit_vals,
      model = sorted_list$code[which(fnames == file)]
    )
  }
) %>%
  bind_rows()
write.csv("summaries/data_pit_st.csv")
data_pit_eta <- read.csv("data_pit_st.csv")

my_palette_0 <- ggsci::pal_lancet()(9)

data_pit_eta %>%
  mutate(model = factor(model, levels = c("NPB", "NEN"))) %>%
  filter(version == "Base") %>%
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
  file.path(fig_path, "fig_9a.eps"),
  width = 3.5,
  height = 3.5
)
# % Reliability diagrams

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
