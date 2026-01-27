# author: Sergio
require(tidyverse)
# auxiliary functions

myposixct <- function(timevar, ...) {
  as.POSIXct(
    case_when(
      nchar(timevar) == 10 ~ paste0(timevar, " 00:00:00"),
      nchar(timevar) == 16 ~ paste0(timevar, ":00"),
      TRUE ~ as.character(timevar)
    ),
    format = "%Y-%m-%d %H:%M:%S",
    ...
  )
}

ar1.trans <- function(x) {
  (exp(x) - 1) / (1 + exp(x))
}


plot.inla.rho.prior <- function(mu = 0, tau = 0.15, n.samples = 1000) {
  # x.test <- seq(-1,1,length.out = 100)
  # plot(x.test, log((x.test+1)/(1-x.test)))

  theta.sample <- rnorm(n.samples, mu, (1 / sqrt(tau)))
  rho.sample <- sapply(theta.sample, ar1.trans)
  # plot(density(theta.sample))
  plot(density(rho.sample))
}


plot.effects <- function(
  inla.model,
  rand.effect,
  trans = \(x) x,
  show.fig = TRUE,
  window = "all",
  ...
) {
  if (rand.effect == "etaderiv" | rand.effect == "eta+deriv") {
    lagged_eta <- inla.model$summary.random$eta.2[
      -nrow(inla.model$summary.random$eta.2),
      c(1, 2, 3, 6, 5, 4, 7, 8)
    ]
    spline_summary <- inla.model$summary.random$eta.1[-1, ] -
      lagged_eta
    # ID fix
    spline_summary$ID <- inla.model$summary.random$eta.1$ID[-1]

    if (rand.effect == "eta+deriv") {
      spline_summary[, 2:ncol(spline_summary)] <- spline_summary[,
        2:ncol(spline_summary)
      ] +
        inla.model$summary.random$eta[-1, 2:ncol(spline_summary)]
    }
  } else {
    # Extract the summary for the specified random effect
    spline_summary <- inla.model$summary.random[[rand.effect]]
  }

  if (window != "all" & is.numeric(window)) {
    spline_summary <- spline_summary %>% tail(trunc(window))
  }

  # Create a data frame for ggplot
  plot_data <- data.frame(
    group = spline_summary$ID, # The grouped values for the random effect
    mean = spline_summary$mean, # Posterior mean
    lower = spline_summary$`0.025quant`, # 2.5% quantile (lower credible interval)
    upper = spline_summary$`0.975quant` # 97.5% quantile (upper credible interval)
  ) %>%
    {
      if (rand.effect %in% c("hour", "month")) {
        mutate(., across(group, as.numeric))
      } else {
        .
      }
    } %>%
    mutate(across(mean:upper, trans))

  # Plot using ggplot2
  p1 <- ggplot(plot_data, aes(x = group, ...)) +
    geom_line(aes(y = mean), color = "blue", lwd = 1) + # Plot mean
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + # Plot credible intervals
    labs(
      title = paste("Estimated effect for", rand.effect),
      x = paste(rand.effect, "(Binned Variable)"),
      y = "Estimated Effect"
    ) +
    theme_minimal()
  if (show.fig) {
    print(p1)
  }
  invisible(list(data = plot_data, fig = p1))
}

plot.hyper.dens <- function(
  inla_model,
  facet.cols = 2,
  logx = FALSE,
  ...
) {
  # Extract the hyperparameter densities into a data frame for plotting
  # browser()
  hyperpar_densities <- do.call(
    rbind,
    lapply(names(inla_model$marginals.hyperpar), function(param) {
      data.frame(
        x = inla_model$marginals.hyperpar[[param]][, "x"],
        y = inla_model$marginals.hyperpar[[param]][, "y"],
        parameter = param
      )
    })
    # mapply(function(param, logtrans) {
    #   data.frame(
    #     x = inla_model$marginals.hyperpar[[param]][, "x"],
    #     y = inla_model$marginals.hyperpar[[param]][, "y"],
    #     parameter = param
    #   ) #%>%
    #     # mutate(logx = ifelse(logtrans, log(x), x))
    # },
    # names(inla_model$marginals.hyperpar),
    # logx,
    # SIMPLIFY = FALSE
    # )
  )
  # browser()
  # hyperpar_densities <- hyperpar_densities %>%
  #   mutate(
  #     x = ifelse(logx, log(x), x)
  #   )

  # Plot the densities using ggplot2
  p.dens <- ggplot(hyperpar_densities, aes(x = x, y = y)) +
    geom_line() +
    facet_wrap(~parameter, scales = "free", ncol = facet.cols) +
    {
      if (logx) {
        scale_x_log10() # Apply log scale to the x-axis
      } else {
        NULL
      }
    } +
    labs(
      title = "Density of Hyperparameters",
      x = "Value",
      y = "Density"
    ) +
    theme_minimal()

  print(p.dens)
}


simulation.plots.inla <- function(
  inla.model,
  data,
  response = "actuals.cf",
  t1 = "2013-06-20 00:01:00 PST",
  quSeq = c(0.025, 0.5, 0.975),
  n.sim.plot = c(3, 100),
  sample.df = NULL,
  family = "beta",
  resp.lab = "Wind Generation",
  ...
) {
  # find points to forecast
  # fcst_points <- which(is.na(train.data[[response]]))
  fcst_dates <- seq(
    from = t1 %>% as.POSIXct(tz = "PST"),
    to = t1 %>% as.POSIXct(tz = "PST") + 24 * 60 * 60,
    by = "hour"
  )
  fcst_points <- which(data$time %in% fcst_dates)
  h <- length(fcst_points)
  nsamp <- 10000

  # Calculate post.pred.samples if sample.df is not provided
  if (is.null(sample.df)) {
    windpow.samples <- inla.posterior.sample(
      n = nsamp,
      result = inla.model,
      selection = list(Predictor = fcst_points),
      num.threads = paste0(mc, ":4")
    )

    precision.samples <- inla.hyperpar.sample(n = nsamp, result = inla.model)
    # Recover errors precision

    if (response == "Y") {
      # reconstruct precision from linear model
      prec.link.npar <- length(inla.model$.args$data$Y) - 2 # excluding response and scale

      scale = inla.model$.args$data$Y$X1 # recover s from data

      phi.samples <- precision.samples[, 1:prec.link.npar] %*% # recover coefficients for precision model
        inla.model$.args$data$Y[2 + 1:prec.link.npar] %>%
        as.data.frame() %>%
        as.matrix() %>%
        t()
    } else {
      # simple precision samples
      phi.samples <- precision.samples[, 1]
    }

    # case logit link
    if (family == "beta") {
      # apply inv.link to samples
      linpredictor.samples <- sapply(windpow.samples, function(x) {
        expit(x$latent)
      }) %>%
        t()

      # recover alpha beta
      shape1.samp <- linpredictor.samples * phi.samples
      shape2.samp <- -linpredictor.samples * phi.samples + phi.samples

      post.pred.samples <- sapply(1:nsamp, \(x) {
        rbeta(h, shape1.samp[x, ], shape2.samp[x, ])
      }) %>%
        t()
    } else {
      # case gaussian, no transformation
      linpredictor.samples <- sapply(windpow.samples, function(x) {
        (x$latent)
      }) %>%
        t()
      # add gaussian noise
      post.pred.samples <- sapply(
        1:nsamp,
        \(x) {
          rnorm(
            h,
            linpredictor.samples[x, ],
            1 / sqrt(phi.samples[x])
          )
        }
      ) %>%
        t()
    }
  } else {
    post.pred.samples <- sample.df
  }

  # Quantiles
  post.pred.quant <- apply(
    post.pred.samples,
    2,
    quantile,
    probs = quSeq,
    na.rm = TRUE
  ) %>%
    t()
  colnames(post.pred.quant) <- paste0("post_pred_q", quSeq)

  # Recover data.frame structure for two linear predictor models
  if (response == "Y") {
    response = attr(inla.model$.args$data$Y, "names.ori")$y # update response
  }

  # Posterior predictive quantiles plot
  plot.q <- data %>%
    filter(date <= "2013-06-21", date >= "2013-06-20") %>%
    cbind(post.pred.quant[(h - 23):h, ]) %>%
    pivot_longer(cols = c(all_of(response), matches("post_pred_q"))) %>%
    ggplot(aes(time, value, col = name)) +
    geom_line() +
    coord_cartesian(...) +
    theme_bw() +
    xlab("date") +
    ylab(paste0(resp.lab, " % of Capacity")) +
    labs(col = "") +
    scale_color_manual(values = c("darkred", rep("gray", length(quSeq)))) +
    theme(legend.position = c(.2, .2))

  # Simulations plot
  select.samples <- post.pred.samples[sample(1:nsamp, n.sim.plot[1]), ] %>% t()
  colnames(select.samples) <- paste0("sim", 1:n.sim.plot[1])
  p.sim.small <- data %>%
    filter(date <= "2013-06-21", date >= "2013-06-20") %>%
    cbind(select.samples[(h - 23):h, ]) %>%
    pivot_longer(cols = c(all_of(response), matches("sim"))) %>%
    mutate(
      name = factor(name, levels = rev(c(response, colnames(select.samples))))
    ) %>%
    ggplot(aes(time, value, col = name)) +
    geom_line() +
    coord_cartesian(...) +
    theme_bw() +
    xlab("date") +
    ylab(paste0(resp.lab, " % of Capacity")) +
    labs(col = "") +
    scale_color_manual(values = rev(c("darkred", rep("gray", n.sim.plot[1])))) +
    theme(legend.position = c(.2, .2))

  select.samples <- post.pred.samples[sample(1:nsamp, n.sim.plot[2]), ] %>% t()
  colnames(select.samples) <- paste0("sim", 1:n.sim.plot[2])
  p.sim.large <- data %>%
    filter(date <= "2013-06-21", date >= "2013-06-20") %>%
    cbind(select.samples[(h - 23):h, ]) %>%
    pivot_longer(cols = c(all_of(response), matches("sim"))) %>%
    mutate(
      name = factor(name, levels = rev(c(response, colnames(select.samples))))
    ) %>%
    ggplot(aes(time, value, col = name)) +
    geom_line() +
    coord_cartesian(...) +
    theme_bw() +
    xlab("date") +
    ylab(paste0(resp.lab, " % of Capacity")) +
    labs(col = "") +
    scale_color_manual(values = rev(c("darkred", rep("gray", n.sim.plot[2])))) +
    theme(legend.position = "none")

  print(plot.q)
  print(p.sim.small)
  print(p.sim.large)
  invisible(post.pred.samples)
}


# plot_actuals_model <- function(
#     data,
#     samples,
#     response = "actuals.cf",
#     t1, # initial time
#     h = 24,
#     resp.lab, # response of the model
#     plot.type = "sim", # simulations or quantiles
#     quSeq = c(0.025,0.5,0.975), # quantiles
#     n.sim.plot = 10, # simulations
#     show.fig = TRUE,
#     legend.opt = "forecast", #
#     ...
# ){
#
#   if (plot.type == "sim"){
#     # filter only a few samples
#     extra.cols <- samples[sample(1:nrow(samples), n.sim.plot), ] %>% t()
#     colnames(extra.cols) <- paste0(plot.type, 1:n.sim.plot)  #add names
#   }
#   if (plot.type == "quant"){
#     # Quantiles
#     extra.cols <- apply(samples, 2, quantile, probs = quSeq, na.rm = TRUE) %>% t()
#     colnames(extra.cols) <- paste0(plot.type, quSeq)
#
#   }
#   plot.data <- data %>%
#     filter(
#       time < t1 + h *60*60,
#       time >= t1) %>%
#     left_join(
#       extra.cols %>% as.data.frame() %>% mutate(time = seq(t1, t1+(nrow(extra.cols)-1)*3600, by ="hour")),
#       by = "time"
#     )
#     # cbind(extra.cols[(h-23):h,])
#
#   model_w_actuals <- model_w_actuals <- (
#     plot.data %>%
#       pivot_longer(cols = c(all_of(response), "forecast.cf", matches(plot.type))) %>%
#       mutate(name = factor(name, levels = rev(c(response, "forecast.cf", colnames(extra.cols))))) %>%
#       {if (response != "actuals.cf") filter(., name != "forecast.cf") else .}
#   ) %>%
#     ggplot(aes(time, value, col = name)) +
#     geom_line() +
#     coord_cartesian(...) +
#     theme_bw() + xlab("date") + ylab(paste0(resp.lab," % of Capacity")) + labs(col = "") +
#     {if (legend.opt=="forecast") {scale_color_manual(
#       values = c("darkred", "darkblue", rep("gray", ncol(extra.cols))),
#         # values = rev(c("darkred", "darkblue", rep("gray", ncol(extra.cols)))),
#         breaks = c("actuals.cf", "forecast.cf") # Exclude names that match plot.type
#       )} else scale_color_manual(
#         values = (c("darkred", rep("gray", ncol(extra.cols)))),
#         breaks = c("err.cf")
#         )}+
#     theme(
#       # legend.position = "bottom",
#       legend.position = "inside",
#       legend.position.inside = c(.3,.2)
#       )
#
#   if(show.fig) print(model_w_actuals)
#   invisible(list(data = plot.data, plot = model_w_actuals))
# }

plot_actuals_model <- function(
  data,
  samples,
  response = "actuals.cf",
  t1, # initial time
  h = 24,
  resp.lab, # response of the model
  plot.type = "sim", # simulations or quantiles
  quSeq = c(0.025, 0.5, 0.975), # quantiles
  n.sim.plot = 10, # simulations
  selection_method = "random",
  show.fig = TRUE,
  legend.opt = "forecast", #
  clipping = FALSE,
  ...
) {
  # recover data if response is inla.mdata

  if (selection_method == "random") {
    # pick random samples
    selection <- sample(1:nrow(samples), n.sim.plot)
  }
  if (selection_method == "rank") {
    # mean power 24 hours
    mpow <- samples %>% apply(., 1, mean)
    mpowrank <- rank(mpow)
    cutoffs <- seq(1, nrow(sampleS), length.out = n.sim.plot + 1) %>%
      round(., digits = 0)
    selection <- mapply(
      \(left, right) {
        # random sample
        # pos <- sample(left:right,size = 1)
        # mid rank
        pos <- round((left + right) / 2, 0)
        # pos
        which(mpowrank == pos)
      },
      left = cutoffs[-(n.sim.plot + 1)],
      right = cutoffs[-1],
      SIMPLIFY = TRUE
    )
  }

  # browser()
  if (plot.type == "sim") {
    # filter only a few samples
    extra.cols <- samples[selection, ] %>% t()
    colnames(extra.cols) <- paste0(plot.type, 1:n.sim.plot) #add names
  }
  if (plot.type == "quant") {
    # Quantiles
    extra.cols <- apply(samples, 2, quantile, probs = quSeq, na.rm = TRUE) %>%
      t()
    colnames(extra.cols) <- paste0(plot.type, quSeq)
  }
  plot.data <- data %>%
    filter(
      time < t1 + h * 60 * 60,
      time >= t1
    ) %>%
    left_join(
      extra.cols %>%
        as.data.frame() %>%
        mutate(time = seq(t1, t1 + (nrow(extra.cols) - 1) * 3600, by = "hour")),
      by = "time"
    ) %>%
    {
      if (response %in% c("err.cf", "Y_err.cf")) {
        # add error samples on top of forecast
        mutate(
          .data = .,
          across(
            c(
              #err.cf,
              all_of(colnames(extra.cols))
            ),
            \(x) forecast.cf_orig + x
          )
        )
      } else {
        .
      }
    } %>%
    {
      if (response %in% c("err", "Y_err")) {
        # add error samples on top of forecast
        mutate(
          .data = .,
          across(
            c(
              #err.cf,
              all_of(colnames(extra.cols))
            ),
            \(x) forecast + x
          )
        )
      } else {
        .
      }
    }

  # cbind(extra.cols[(h-23):h,])
  # browser()
  fcst.name.plot <- case_when(
    response %in% c("actuals", "err", "Y_actuals", "Y_err") ~ "forecast",
    TRUE ~ "forecast.cf_orig"
  )
  obs.name.plot <- case_when(
    response %in% c("actuals", "err", "Y_actuals", "Y_err") ~ "actuals",
    TRUE ~ "actuals.cf"
  )
  updated_response <- if (grepl("Y_", response)) {
    substr(response, 3, nchar(response))
  } else {
    response
  }
  model_w_actuals <- (
    plot.data %>%
      pivot_longer(
        cols = c(
          all_of(updated_response),
          all_of(c(fcst.name.plot, obs.name.plot)),
          matches(plot.type)
        )
      ) %>%
      mutate(
        name = factor(
          name,
          levels = rev(c(
            updated_response,
            fcst.name.plot,
            obs.name.plot,
            colnames(extra.cols)
          )) %>%
            unique()
        )
      ) %>%
      {
        if (grepl("err", response)) filter(., !grepl("err", name)) else .
      } %>% # remove observed error
      # {if (response %in% c("err.cf","Y_err.cf")) filter(., name != "err.cf") else .} %>%  # remove error
      {
        if (response %in% c("actuals", "Y_actuals")) {
          filter(., name != "actuals.cf")
        } else {
          .
        }
      } # remove normalised
  ) %>%
    {
      if (clipping) mutate(., value = pmin(1, pmax(0, value))) else .
    } %>%
    ggplot(aes(time, value, col = name)) +
    geom_line() +
    coord_cartesian(...) +
    theme_bw() +
    # xlab("date") +
    ylab(paste0(resp.lab, " % of Capacity")) +
    labs(col = "") +
    # scale_color_manual(values = rep("lightblue", ncol(extra.cols)+2)) +
    {
      if (legend.opt == "normalised") {
        scale_color_manual(
          values = c("darkred", "darkblue", rep("lightblue", ncol(extra.cols))),
          breaks = c("actuals.cf", "forecast.cf_orig"), # Exclude names that match plot.type
          labels = c("observed", "forecast")
        )
      } else {
        if (legend.opt == "error") {
          scale_color_manual(
            values = (c(
              "darkred",
              "darkblue",
              rep("lightblue", ncol(extra.cols))
            )),
            breaks = c("actuals.cf", "forecast.cf_orig"),
            labels = c("observed", "forecast")
          )
        } else {
          scale_color_manual(
            values = c(
              "darkred",
              "darkblue",
              rep("lightblue", ncol(extra.cols))
            ),
            breaks = c("actuals", "forecast"), # Exclude names that match plot.type
            labels = c("observed", "forecast")
          )
          # scale_color_manual(
          #   values = (c("darkred", rep("lightblue", ncol(extra.cols)))),
          #   breaks = c("err.cf"),
          #   labels = c("forecast error")
          # )
        }
      }
    } +
    labs(x = "") +
    scale_x_datetime(date_labels = "%H:%M") +
    theme(
      legend.position = "none",
      # legend.position = "bottom",
      # legend.position = "inside",
      legend.position.inside = c(.2, .8),
      legend.background = element_blank(), # Makes background completely transparent
      legend.box.background = element_rect(fill = NA, color = NA) # No border
    )

  if (show.fig) {
    print(model_w_actuals)
  }
  invisible(list(data = plot.data, plot = model_w_actuals))
}

simulation.plots.inla2 <- function(
  inla.model,
  data,
  response = "actuals.cf",
  t1 = "2013-06-20 00:00:00 PST",
  h = 24,
  quSeq = c(0.025, 0.5, 0.975),
  n.sim.plot = c(3, 100),
  sample.df = NULL,
  nsamp = 10000,
  family = "beta",
  resp.lab = "Wind Generation",
  show.fig = TRUE,
  save.fig = FALSE,
  path = NULL,
  run.name = "default",
  skip.plots = FALSE,
  inla_seed = 1,
  ...
) {
  # browser()
  # Conditionally convert to POSIXct only if t1 is character
  t1 <- t1 %>%
    {
      if (is.character(.)) {
        as.POSIXct(., format = "%Y-%m-%d %H:%M:%S", tz = "PST")
      } else {
        .
      }
    }
  if (grepl("Y_", response) & !is.na(inla.model$.args$data$time[1])) {
    pos_shift <- length(inla.model$.args$data$time) / 2
  } else {
    pos_shift <- 0
  }
  # find points to forecast
  fcst_dates <- seq(
    from = t1 + 0 * 60 * 60,
    to = t1 + h * 60 * 60,
    by = "hour"
  )
  fcst_points <- which(inla.model$.args$data$time %in% fcst_dates) + pos_shift
  # h <- length(fcst_points)

  # Calculate post.pred.samples if sample.df is not provided
  if (is.null(sample.df)) {
    windpow.samples <- inla.posterior.sample(
      n = nsamp,
      result = inla.model,
      selection = list(Predictor = fcst_points),
      # num.threads = paste0(mc, ":4")
      # use.improved.mean = FALSE,
      # verbose = TRUE,
      num.threads = "1:1",
      seed = inla_seed
    )

    precision.samples <- inla.hyperpar.sample(n = nsamp, result = inla.model)
    # Recover errors precision

    hypers <- inla.model$summary.hyperpar %>% rownames()

    # number of ar terms
    ar.order <- max(1, sum(str_detect(hypers, "^PACF[0-9]+ for t$")))

    # find out which ar model was used
    ar.type <- case_when(
      any(grepl("Rho for t", hypers)) ~ "AR1",
      any(grepl("PACF1 for t", hypers)) ~ paste0("AR", ar.order),
      TRUE ~ "none"
    )
    # browser()

    if (ar.type == "AR1") {
      # get ar par samples
      sigma_ar <- precision.samples[, "Precision for t"]
      rho <- precision.samples[, "Rho for t"]
    } else {
      if (ar.type != "none") {
        # get ar par samples
        sigma_ar <- precision.samples[, "Precision for t"]
        position <- which(grepl(
          pattern = "^PACF[0-9]+ for t$",
          colnames(precision.samples)
        ))
        pacf_vec <- precision.samples[, position]
      }
    }

    if (response == "Y") {
      # for generalised Gaussian
      # reconstruct precision from linear model
      prec.link.npar <- length(inla.model$.args$data$Y) - 2 # excluding response and scale

      scale = inla.model$.args$data$Y$X1[fcst_points] # recover s from data

      # recover coefficients for precision model
      phi.samples <- precision.samples[, 1:prec.link.npar] %*% # nsamp x npar
        (inla.model$.args$data$Y[2 + 1:prec.link.npar] %>%
          as.data.frame() %>%
          slice(fcst_points) %>% # only times to forecast
          as.matrix() %>%
          t()) # n.times x npar

      phi.samples <- scale * exp(phi.samples) # convert to original scale
    } else {
      # simple precision samples
      phi.samples <- precision.samples[, 1]
    }
    # case logit link
    if (family == "beta") {
      # apply inv.link to samples
      linpredictor.samples <- sapply(windpow.samples, function(x) {
        expit(x$latent)
      }) %>%
        t()

      # recover alpha beta
      shape1.samp <- linpredictor.samples * phi.samples
      shape2.samp <- -linpredictor.samples * phi.samples + phi.samples

      post.pred.samples <- sapply(
        1:nsamp,
        \(x) {
          rbeta(ncol(linpredictor.samples), shape1.samp[x, ], shape2.samp[x, ])
        }
      ) %>%
        t()
      # parametric_quants <- sapply(
      #   quSeq,
      #   \(prob) qbeta(p = prob,
      #                 shape1 = shape1.samp,
      #                 shape2 = shape2.samp)
      # )
      # colnames(parametric_quants) <- paste0("quant",probs)
    } else {
      # case gaussian, no transformation
      linpredictor.samples <- sapply(
        windpow.samples,
        function(x) {
          if (family %in% c("stochvol", "stochvolt", "gamma", "weibull")) {
            exp(x$latent)
          } else {
            x$latent
          }
        }
      ) %>%
        t()
      # browser()
      if (family %in% c("gamma")) {
        post.pred.samples <- sapply(
          1:nsamp,
          \(x) {
            rgamma(
              ncol(linpredictor.samples),
              phi.samples[x], # alpha
              phi.samples[x] / linpredictor.samples[x, ] # rate
            )
          }
        ) %>%
          t()
      } else if (family %in% c("weibull")) {
        post.pred.samples <- sapply(
          1:nsamp,
          \(x) {
            rweibull(
              ncol(linpredictor.samples),
              phi.samples[x], # alpha
              # 1/linpredictor.samples[x,] # scale 1/lambda
              linpredictor.samples[x, ]^(-1 / phi.samples[x]) # scale lambda^(-1/alpha)
            )
          }
        ) %>%
          t()
      } else {
        post.pred.samples <- sapply(
          1:nsamp,
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
      }
      # add gaussian noise
      # post.pred.samples <- sapply(
      #   1:nsamp,
      #   \(x) rnorm(
      #     ncol(linpredictor.samples),
      #     linpredictor.samples[x,],1/sqrt(phi.samples[x]))
      # ) %>% t()
    }
  } else {
    post.pred.samples <- sample.df
  }

  # Recover data.frame structure for two linear predictor models
  if (response == "Y") {
    # for generalised gaussian
    response = attr(inla.model$.args$data$Y, "names.ori")$y # update response
    # for eta derivative trick
  }
  # browser()

  if (!skip.plots) {
    legend.opt <- case_when(
      grepl("actuals.cf", response) ~ "normalised",
      grepl("err.cf", response) ~ "error",
      grepl("actuals", response) ~ "power",
      grepl("err", response) ~ "err"
    )
    plot.q <- plot_actuals_model(
      data,
      samples = post.pred.samples,
      response = response,
      t1 = t1,
      h = h,
      resp.lab = resp.lab,
      plot.type = "quant",
      quSeq = quSeq,
      # n.sim.plot = n.sim.plot[1],
      show.fig = FALSE,
      legend.opt = legend.opt,
      ...
    )

    p.sim.small <- plot_actuals_model(
      data,
      samples = post.pred.samples,
      response = response,
      t1 = t1,
      h = h,
      resp.lab = resp.lab,
      plot.type = "sim",
      n.sim.plot = n.sim.plot[1],
      show.fig = FALSE,
      legend.opt = legend.opt,
      ...
    )
    p.sim.large <- plot_actuals_model(
      data,
      samples = post.pred.samples,
      response = response,
      t1 = t1,
      h = h,
      resp.lab = resp.lab,
      plot.type = "sim",
      n.sim.plot = n.sim.plot[2],
      show.fig = FALSE,
      legend.opt = legend.opt,
      ...
    )

    if (show.fig) {
      print(plot.q$plot)
      print(p.sim.small$plot)
      print(p.sim.large$plot)
    }
    if (save.fig) {
      # Create the directory if it doesn't exist
      if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
      }
      ggsave(
        file.path(path, paste0(run.name, "_quantiles.eps")),
        plot.q$plot,
        width = 3.5,
        height = 3
      )
      ggsave(
        file.path(path, paste0(run.name, "_sim.small.eps")),
        p.sim.small$plot,
        width = 3.5,
        height = 3
      )
      ggsave(
        file.path(path, paste0(run.name, "_sim.large.eps")),
        p.sim.large$plot,
        width = 3.5,
        height = 3
      )
    }
  } else {
    plot.q <- p.sim.small <- p.sim.large <- list()
    plot.q$data <- NULL
    plot.q$plot <- NULL
    p.sim.small$plot <- NULL
    p.sim.large <- NULL
  }
  invisible(
    list(
      samples = post.pred.samples,
      quantiles = plot.q$data,
      plot.q = plot.q$plot,
      p.sim.small = p.sim.small$plot,
      p.sim.large = p.sim.large$plot
    )
  )
}

get_inla_formula <- function(inla_model) {
  as.formula(
    inla_model$.args$formula %>% deparse(backtick = T) %>% paste(collapse = "")
  )
}

# model posterior estimation for a given date
# fit_a_date <- function(
#     timet = "2013-06-20 00:01:00 PST",
#     h = 24,
#     dat.power, # should have time, actuals.cf, ws, month, hour,
#     response,
#     timezone = "PST",
#     inla.object = NULL,
#     restart = FALSE,
#     ...){
#
#   # dates of interest
#   time0 <- min(dat.power$time) #, initial time,
#   timet <- as.POSIXct(timet, tz=timezone) # last time included in training
#   timeth <- timet + as.difftime(h, units = "hours") # last time to predict
#
#   # prepare data
#   model_data <- dat.power %>%
#     # mutate(actuals.cf = na.approx(actuals.cf)) %>%
#     # simple time index
#     # mutate(t = as.numeric(difftime(time, first(time), units = "hours"))) %>%
#     filter(
#       time >= time0,
#       time <=timeth) %>%
#     # mutate(
#     #   month= month(time) %>% factor(),
#     #   hour = hour(time) %>% factor()
#     # ) %>%
#     # set target response to NA to be able to forecast it
#     mutate(across(c(actuals.cf,actuals,err.cf), \(x) ifelse(time < timet, x , NA)))
#     # mutate(actuals.cf = ifelse(time < timet, actuals.cf,NA)) %>%
#     # mutate(err.cf = ifelse(time < timet, err.cf, NA))
#
#   # adjustments in case of ggaussian family
#   if(inla.object$.args$family == "ggaussian"){
#     s = rep(1,nrow(model_data))
#     interc = rep(1,nrow(model_data))
#     Y = with(
#       model_data,
#       inla.mdata(err.cf, s, interc, forecast.cf, forecast.cf^2))
#   }
#
#   # data
#   if (inla.object$.args$family != "ggaussian") {
#     day.data <- model_data
#   } else {
#     day.data <- list(
#       Y = Y,
#       ws.w_group = model_data$ws.w_group,
#       t = model_data$t,
#       month = model_data$month,
#       hour = model_data$hour,
#       intercept = interc,
#       forecast.cf = model_data$forecast.cf,
#       time = model_data$time,
#       s = s
#     )
#   }
#
#
#   # prepare model run
#   new_inla_model <- inla(
#     formula = get_inla_formula(inla.object), # model formula
#     data = day.data,
#     family = inla.object$.args$family,
#     control.mode = list(theta = inla.object$mode$theta, ...), # do not restart opt
#     control.family = inla.object$.args$control.family,
#     control.compute = inla.object$.args$control.compute,
#     control.predictor = list(compute = inla.object$.args$control.predictor$compute,
#                              link = 1),
#     # verbose = TRUE,
#     ...
#   )
#   # return INLA object
#   new_inla_model
# }
fit_a_date <- function(
  timet = "2013-06-20 00:01:00 PST",
  h = 24,
  dat.power, # should have time, actuals.cf, ws, month, hour,
  response,
  timezone = "PST",
  inla.object = NULL,
  restart = FALSE,
  ...
) {
  # dates of interest
  time0 <- min(dat.power$time) #, initial time,
  timet <- as.POSIXct(timet, tz = timezone) # last time included in training
  timeth <- timet + as.difftime(h, units = "hours") # last time to predict
  # browser()
  # prepare data
  model_data <- dat.power %>%
    # mutate(actuals.cf = na.approx(actuals.cf)) %>%
    # simple time index
    # mutate(t = as.numeric(difftime(time, first(time), units = "hours"))) %>%
    filter(
      time >= time0,
      time <= timeth
    ) %>%
    # mutate(
    #   month= month(time) %>% factor(),
    #   hour = hour(time) %>% factor()
    # ) %>%
    # set target response to NA to be able to forecast it
    mutate(across(c(actuals.cf, actuals, err.cf, err), \(x) {
      ifelse(time < timet, x, NA)
    }))
  # mutate(actuals.cf = ifelse(time < timet, actuals.cf,NA)) %>%
  # mutate(err.cf = ifelse(time < timet, err.cf, NA))

  # exclude zeroes
  if (tail(inla.object$.args$family, 1) %in% c("weibull", "gamma")) {
    model_data <- model_data %>%
      filter(is.na(actuals) | actuals > 0)
  }

  # adjustments in case of etaderiv
  if (length(inla.object$.args$family) > 1) {
    n <- nrow(model_data)
    updated_response <- gsub("Y_", "", response)
    # change data
    day.data <- with(
      model_data,
      list(
        # Dynamically set the name using updated_response
        intercept = c(rep(1, n), rep(NA, n)), # 1st likelihood intercept, nothing for 2nd as it's stored in eta
        # wind = c(ws.w, rep(NA, n)), # x data goes in 1st likelihood
        ws.w_group = c(ws.w_group, rep(NA, n)),
        fcst_group = c(fcst_group, rep(NA, n)),
        fd_group = c(fd_group, rep(NA, n)),
        t = c(t, rep(NA, n)),
        time = c(rep(NA, n), time),
        month = c(month, rep(NA, n)),
        hour = c(hour, rep(NA, n)),
        eta = c(1:n, 1:n), # eta indices
        w = c(rep(-1, n), rep(1, n)), # weights for eta effect: -1 to copy lin.predictor, 1 to be part of likelihood
        eta.1 = c(rep(NA, n), 1:n), # indices for positive part of derivative
        eta.2 = c(rep(NA, n), NA, 1:(n - 1)), # shifted indices for negative part of derivative
        w2 = c(rep(0, n), rep(-1, n))
      ) # weights for negative part of derivative
    )
    # adding properly named response
    day.data[[response]] <- cbind(
      c(rep(0, n), rep(NA, n)), # fake zeros (1st likelihood)
      c(rep(NA, n), model_data[[updated_response]])
    )
  } else {
    # adjustments in case of ggaussian family
    if (inla.object$.args$family == "ggaussian") {
      s = rep(1, nrow(model_data))
      interc = rep(1, nrow(model_data))
      Y = with(
        model_data,
        inla.mdata(err.cf, s, interc, forecast.cf, forecast.cf^2)
      )
    }
    # data
    if (inla.object$.args$family != "ggaussian") {
      day.data <- model_data
    } else {
      day.data <- list(
        Y = Y,
        ws.w_group = model_data$ws.w_group,
        t = model_data$t,
        month = model_data$month,
        hour = model_data$hour,
        intercept = interc,
        forecast.cf = model_data$forecast.cf,
        time = model_data$time,
        s = s
      )
    }
  }
  # browser()
  # prepare model run
  new_inla_model <- inla(
    formula = get_inla_formula(inla.object), # model formula
    data = day.data,
    family = inla.object$.args$family,
    control.mode = list(theta = inla.object$mode$theta, restart = restart), # do not restart opt
    control.family = inla.object$.args$control.family,
    control.compute = inla.object$.args$control.compute,
    control.predictor = list(
      compute = inla.object$.args$control.predictor$compute,
      link = 1
    ),
    # verbose = TRUE,
    ...
  )
  # return INLA object
  new_inla_model
}

all.visualistions <- function(
  inla_model,
  t1,
  data,
  effects = TRUE,
  samples = TRUE,
  hyper = TRUE,
  ...
) {
  if (effects) {
    # Effects
    effects.list <- names(inla_model$summary.random)
    excluded <- c("t", "eta", "eta.1", "eta.2")

    lapply(
      effects.list[!effects.list %in% excluded],
      \(effect) plot.effects(inla_model, effect)
    )
  }
  # plot.effects(inla_model,"ws.w_group")
  # plot.effects(inla_model,"fcst_group")
  # plot.effects(mod.test,"month", trans = \(x) x, group =1)
  # plot.effects(mod.test,"hour", trans = \(x) x, group =1)

  if (samples) {
    # Samples
    sample.test <- simulation.plots.inla2(
      inla.model = inla_model,
      data = data,
      response = inla_model$.args$formula[[2]] %>% as.character(),
      t1 = t1,
      quSeq = c(0.025, 0.5, 0.975),
      family = tail(inla_model$.args$family, 1),
      resp.lab = "Wind Power",
      # ylim = c(0,1.02),
      nsamp = 1000,
      # sample.df = sample.test$samples,
      # show.fig = FALSE,
      # save.fig = TRUE,
      # path = file.path("~/Documents/proj2","fig"),
      # run.name = paste0("r_act_f_beta_t",t1)
      # legend.position = "bottom"
      ...
    )
  }
  # undebug(simulation.plots.inla2)

  if (hyper) {
    # Densities
    plot.hyper.dens(inla_model, ...)
  }
}


# functions_to_scenarios <- function(data, func_list, ...) {
#   # Check that func_list is a list of functions
#   if (!all(sapply(func_list, is.function))) {
#     stop("func_list should be a list of functions.")
#   }
#
#   # # Apply each function to each column and store the results in a list of data frames
#   # results <- lapply(func_list, function(f) {
#   #   # Apply function f to each column (time point)
#   #   apply(data, 2, f, ...)
#   # }, ...)
#
#   # Apply each function to each column and store the results in a list of data frames
#   results <- lapply(func_list, function(f) {
#     # Check if the function can take additional arguments by checking its formal arguments
#     func_args <- names(formals(f))
#     if ("..." %in% func_args || all(names(list(...)) %in% func_args)) {
#       # Pass extra arguments only if the function accepts them
#       apply(data, 2, function(x) f(x, ...))
#     } else {
#       # If it doesn’t accept additional arguments, call it without ...
#       apply(data, 2, f)
#     }
#   })
#
#   # Combine results into a single data.frame, each function result as a new column
#   result_df <- as.data.frame(do.call(cbind, results))
#
#   # Rename columns with function names for clarity
#   # names(result_df) <- sapply(func_list, function(f) deparse(substitute(f)))
#
#   # Set row names as time points (1 to 24) and reset row names for tidy output
#   # rownames(result_df) <- 1:ncol(data)
#   result_df <- tibble::rownames_to_column(result_df, var = "time")
#
#   return(result_df)
# }

functions_to_scenarios <- function(data, func_list, ...) {
  # Check that func_list is a list of functions
  if (!all(sapply(func_list, is.function))) {
    stop("func_list should be a list of functions.")
  }

  # Apply each function to each column and store the results in a list
  results <- lapply(func_list, function(f) {
    # Check if the function can take additional arguments by checking its formal arguments
    func_args <- names(formals(f))
    if ("..." %in% func_args || all(names(list(...)) %in% func_args)) {
      # Apply function with extra arguments and check if result has multiple values
      apply(data, 2, function(x) f(x, ...))
    } else {
      # Apply function without extra arguments
      apply(data, 2, f)
    }
  })

  # Expand the results list, handling functions that return multiple values
  expanded_results <- lapply(seq_along(results), function(i) {
    res <- results[[i]]
    if (is.matrix(res)) {
      # If a matrix (e.g., multiple quantiles), convert to data.frame with specific column names
      as.data.frame(t(res))
    } else {
      # Otherwise, just convert to a data.frame
      data.frame(res)
    }
  })

  # Combine all expanded results into a single data frame, naming columns appropriately
  result_df <- do.call(cbind, expanded_results)

  # Set column names for each function's output (e.g., mean, sd, quant1, quant2, quant3)
  func_names <- names(func_list)
  col_names <- unlist(lapply(seq_along(func_list), function(i) {
    if (is.matrix(results[[i]])) {
      # paste(func_names[i], seq_len(nrow(results[[i]])), sep = "_")
      paste(func_names[i], rownames(results[[i]]), sep = "_")
    } else {
      func_names[i]
    }
  }))

  names(result_df) <- col_names
  # rownames(result_df) <- 1:ncol(data)

  # Reset row names for tidy output, adding time points
  result_df <- tibble::rownames_to_column(result_df, var = "time")

  return(result_df)
}


scores_day_model <- function(
  sample.path = file.path("~/Documents/proj2", "sample"),
  model.name,
  t1,
  response
) {
  scen.tbl <- file.path(
    sample.path,
    paste0(model.name, date(t1), ".csv")
  ) %>%
    fread(.) %>%
    as.data.frame()
  # compare with actuals - forecast error
  fcst_dates <- seq(
    from = t1 + 0 * 60 * 60,
    to = t1 + 23 * 60 * 60,
    by = "hour"
  )

  subsample <- data.scaled %>%
    filter(time %in% fcst_dates)

  # orig.obs <- subsample %>%
  #   select(time,all_of(c("","forecast.cf")))

  # Define functions to apply (e.g., mean, sd, quantile)
  function_list <- list(mean = mean, sd = sd, quantile = quantile)

  # Call the function with desired quantiles
  # result <- functions_to_scenarios(sim_data, function_list, probs = c(0.05, 0.5, 0.95))
  result <- functions_to_scenarios(
    scen.tbl,
    function_list,
    probs = c(0.005, 0.025, 0.5, 0.975)
  ) %>%
    mutate(
      # Remove "t_" prefix and convert to POSIXct, handling NA text entries
      time = if_else(
        str_detect(time, "^NA$"), # Check if the value is "NA"
        NA_POSIXct_, # Keep it as NA if it's "NA"
        as.POSIXct(str_remove(time, "^t_"), format = "%Y-%m-%d_%H", tz = "PST")
      )
    )

  checks <- subsample %>%
    left_join(result, by = "time") %>%
    # mutate(coverage = between(.data[[response]],`quantile_2.5%`,`quantile_97.5%`))
    summarise(
      rmse = sqrt(mean((.data[[response]] - mean)^2)),
      coverage95 = mean(between(
        .data[[response]],
        `quantile_2.5%`,
        `quantile_97.5%`
      )),
      coverage005 = mean(.data[[response]] < `quantile_0.5%`)
    ) %>%
    mutate(date = t1)
  return(checks)
}


# model.scoring <- function(
#     model.name,
#     time_seq,
#     sample.path = file.path("~/Documents/proj2","sample")
# ){
#   # get response from model code
#   response = if_else(
#     str_detect(model.name, "r_act_f"),
#     "actuals.cf",
#     if_else(str_detect(model.name, "r_err_f"), "err.cf", NA_character_)
#   )
#
#   # get daily scores
#   scores.model <- lapply(
#     time_seq,
#     \(day) scores_day_model(
#       sample.path ,
#       model.name,
#       day,
#       response)) %>%
#     do.call(bind_rows,.)
#
#   # get mean scores
#   mean.scores <- scores.model %>%
#     # do.call(bind_rows,.) %>%
#     summarise(across(where(is.numeric), \(x) mean(x,na.rm = TRUE))) %>%
#     mutate(
#       model = model.name,
#       response = response
#     )
#
#   return(list(day = mean.scores, hour = scores.model))
#
# }

### new hourly stats
stats_hour_model <- function(
  sample.path = file.path("~/Documents/proj2", "sample"),
  model.name,
  t1,
  response,
  probs = c(0.005, 0.025, 0.5, 0.975),
  compressed = FALSE,
  clipping = FALSE,
  power_data = data.scaled
) {
  if (compressed) {
    fname <- paste0(model.name, date(t1), ".csv.gz")
  } else {
    fname <- paste0(model.name, date(t1), ".csv")
  }

  # browser()
  scen.tbl <- file.path(
    sample.path,
    fname
  ) %>%
    fread(.) %>%
    as.data.frame()

  # compare with actuals - forecast error
  fcst_dates <- seq(
    from = t1 + 0 * 60 * 60,
    to = t1 + 23 * 60 * 60,
    by = "hour"
  )

  subsample <- power_data %>%
    filter(time %in% fcst_dates)

  # browser()
  additive_factor <- case_when(
    grepl("err.cf$", response) ~ subsample[["forecast.cf_orig"]],
    grepl("err$", response) ~ subsample[["forecast"]],
    TRUE ~ 0,
  )

  mult_factor <- case_when(
    grepl("actuals$", response) ~ subsample[["capacity"]],
    grepl("err$", response) ~ subsample[["capacity"]],
    TRUE ~ 1,
  )
  # rescaling
  scen.tbl <- apply(
    scen.tbl[, 1:24],
    1,
    \(row) {
      rescaled_sample <- (row + additive_factor) / mult_factor
      if (clipping) {
        rescaled_sample <- pmin(1, pmax(0, rescaled_sample))
      }
      rescaled_sample
    }
  ) %>%
    t()
  observed <- "actuals.cf"
  # orig.obs <- subsample %>%
  #   select(time,all_of(c("","forecast.cf")))

  # Define functions to apply (e.g., mean, sd, quantile)
  function_list <- list(
    mean = mean,
    sd = sd,
    quantile = quantile
  )

  # Call the function with desired quantiles
  # result <- functions_to_scenarios(sim_data, function_list, probs = c(0.05, 0.5, 0.95))
  result <- functions_to_scenarios(scen.tbl, function_list, probs = probs) %>%
    mutate(
      # Remove "t_" prefix and convert to POSIXct, handling NA text entries
      # time = if_else(
      #   str_detect(time, "^NA$"), # Check if the value is "NA"
      #   NA_POSIXct_,              # Keep it as NA if it's "NA"
      #   as.POSIXct(str_remove(time, "^t_"), format = "%Y-%m-%d_%H", tz = "PST")
      # )
      time = fcst_dates
    )

  crps <- scoringRules::crps_sample(
    y = subsample %>% pull(!!observed), # day-ahead response
    dat = scen.tbl %>%
      as.data.frame() %>%
      select(1:nrow(subsample)) %>% # select 24-hours
      as.matrix() %>%
      t() # transpose so that crps works
  )

  energy.s <- scoringRules::es_sample(
    y = subsample %>% pull(!!observed), # day-ahead response
    dat = scen.tbl %>%
      as.data.frame() %>%
      select(1:nrow(subsample)) %>% # select 24-hours
      as.matrix() %>%
      t() # transpose so that crps works
  )

  # variogram distance
  p.var <- c(0.5, 1, 2)
  t.ahead <- c(1, seq(6, 24, by = 6))
  # variogram <- sapply(
  #   p.var,
  #   \(p) scoringRules::vs_sample(
  #     y = subsample %>%
  #       slice(t.ahead) %>%
  #       pull(!!response), # day-ahead response
  #     dat = scen.tbl %>%
  #       select(all_of(t.ahead)) %>% # select times
  #       as.matrix() %>% t(), # transpose so that crps works
  #     p = p
  #   ))

  # average of the day
  checks <- subsample %>%
    left_join(result, by = "time") %>%
    # mutate(coverage = between(.data[[response]],`quantile_2.5%`,`quantile_97.5%`))
    #   summarise(
    #     rmse = sqrt(mean((.data[[response]] - mean)^2)),
    #     coverage95 = mean(between(.data[[response]],`quantile_2.5%`,`quantile_97.5%`)),
    #     coverage005 = mean(.data[[response]]<`quantile_0.5%`)
    #   ) %>%
    mutate(
      date = t1,
      crps = crps,
      energy = energy.s, # recycling: energy score is just 1 number per day-ahed fcst.
      !!!setNames(
        lapply(p.var, function(p) {
          scoringRules::vs_sample(
            y = subsample %>%
              slice(t.ahead) %>%
              pull(!!observed), # day-ahead response
            dat = scen.tbl %>%
              as.data.frame() %>%
              select(all_of(t.ahead)) %>% # select times
              as.matrix() %>%
              t(), # Transpose for vs_sample
            p = p
          )
        }),
        paste0("variogram_p_", p.var)
      )
    )
  return(checks)
}

model.scoring <- function(
  model.name,
  time_seq,
  sample.path = file.path("~/Documents/proj2", "sample"),
  probs = c(0.005, 0.025, 0.5, 0.975)
) {
  # get response from model code
  response = if_else(
    str_detect(model.name, "r_act_f"),
    "actuals.cf",
    if_else(str_detect(model.name, "r_err_f"), "err.cf", NA_character_)
  )

  # get hourly stats
  hour_stats <- lapply(
    time_seq,
    \(day) {
      stats_hour_model(
        sample.path,
        model.name,
        day,
        response,
        probs
      )
    }
  ) %>%
    do.call(bind_rows, .)

  daily.scores.model <- hour_stats %>%
    group_by(date) %>%
    summarise(
      rmse = sqrt(mean((.data[[response]] - mean)^2)),
      coverage95 = mean(between(
        .data[[response]],
        `quantile_2.5%`,
        `quantile_97.5%`
      )),
      coverage005 = mean(.data[[response]] < `quantile_0.5%`)
    ) %>%
    mutate(
      model = model.name,
      response = response
    )

  hour.scores.model <- hour_stats %>%
    group_by(hour) %>%
    summarise(
      rmse = sqrt(mean((.data[[response]] - mean)^2)),
      coverage95 = mean(between(
        .data[[response]],
        `quantile_2.5%`,
        `quantile_97.5%`
      )),
      coverage005 = mean(.data[[response]] < `quantile_0.5%`)
    ) %>%
    mutate(
      model = model.name,
      response = response
    )

  # get mean scores
  mean.scores <- daily.scores.model %>%
    # do.call(bind_rows,.) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
    mutate(
      model = model.name,
      response = response
    )

  return(list(
    global = mean.scores,
    day = daily.scores.model,
    hour = hour.scores.model
  ))
}

model.scoring.reliability <- function(
  model.name,
  time_seq,
  sample.path = file.path("~/Documents/proj2", "sample"),
  compressed = FALSE,
  parallel = FALSE,
  mc = 1,
  ...
) {
  # browser()
  # get response from model code
  # response = if_else(
  #   str_detect(model.name, "r_act_f|r_actuals-cf_f"),
  #   "actuals.cf",
  #   if_else(str_detect(model.name, "r_err_f"), "err.cf", NA_character_)
  # )
  response <- case_when(
    str_detect(model.name, "r_act_f|r_actuals-cf_f") ~ "actuals.cf",
    TRUE ~ str_extract(model.name, "(?<=r_).*?(?=_f)") %>%
      gsub("-", ".", .)
  )

  probs = seq(0.05, 0.95, by = 0.05) # probs for reliability diagram
  alpha <- 1 - probs # alpha level
  quant_vector <- c(alpha / 2, 1 - alpha / 2) # quantiles for symmetric confidence interval

  # get hourly stats
  # hour_stats <- lapply(
  #   time_seq,
  #   \(day) stats_hour_model(
  #     sample.path,
  #     model.name,
  #     day,
  #     response,
  #     probs = c(quant_vector, 0.005), # adding half a percent for extreme low wind
  #     compressed = compressed
  #   )) %>%
  #   do.call(bind_rows,.)
  # browser()
  if (parallel) {
    hour_stats <- parallel::mclapply(
      time_seq,
      \(day) {
        tryCatch(
          stats_hour_model(
            sample.path,
            model.name,
            day,
            response,
            probs = c(quant_vector, 0.005),
            compressed = compressed,
            ...
          ),
          error = function(e) data.frame() # Return an empty data.frame() if an error occurs
        )
      },
      mc.cores = mc
    ) %>%
      bind_rows()
  } else {
    hour_stats <- lapply(
      time_seq,
      \(day) {
        tryCatch(
          stats_hour_model(
            sample.path,
            model.name,
            day,
            response,
            probs = c(quant_vector, 0.005),
            compressed = compressed,
            ...
          ),
          error = function(e) data.frame() # Return an empty data.frame() if an error occurs
        )
      }
    ) %>%
      bind_rows()
  }

  # browser()
  # quantile related column names
  quant_low <- paste0("quantile_", round(alpha / 2 * 100, 1), "%")
  quant_high <- paste0("quantile_", round((1 - alpha / 2) * 100, 1), "%")

  stats.vec <- c(
    "mean",
    "sd",
    "crps",
    "energy",
    "variogram_p_0.5",
    "variogram_p_1",
    "variogram_p_2"
  )
  observed <- "actuals.cf"
  # calculate all coverages by day
  daily.scores.model <- hour_stats %>%
    group_by(date) %>%
    summarise(
      rmse = sqrt(mean((.data[[observed]] - mean)^2)),
      across(!!stats.vec, mean),
      # crps = mean(crps)
      coverage05 = mean(between(
        .data[[observed]],
        .data[[quant_low[1]]],
        .data[[quant_high[1]]]
      )),
      coverage10 = mean(between(
        .data[[observed]],
        .data[[quant_low[2]]],
        .data[[quant_high[2]]]
      )),
      coverage15 = mean(between(
        .data[[observed]],
        .data[[quant_low[3]]],
        .data[[quant_high[3]]]
      )),
      coverage20 = mean(between(
        .data[[observed]],
        .data[[quant_low[4]]],
        .data[[quant_high[4]]]
      )),
      coverage25 = mean(between(
        .data[[observed]],
        .data[[quant_low[5]]],
        .data[[quant_high[5]]]
      )),
      coverage30 = mean(between(
        .data[[observed]],
        .data[[quant_low[6]]],
        .data[[quant_high[6]]]
      )),
      coverage35 = mean(between(
        .data[[observed]],
        .data[[quant_low[7]]],
        .data[[quant_high[7]]]
      )),
      coverage40 = mean(between(
        .data[[observed]],
        .data[[quant_low[8]]],
        .data[[quant_high[8]]]
      )),
      coverage45 = mean(between(
        .data[[observed]],
        .data[[quant_low[9]]],
        .data[[quant_high[9]]]
      )),
      coverage50 = mean(between(
        .data[[observed]],
        .data[[quant_low[10]]],
        .data[[quant_high[10]]]
      )),
      coverage55 = mean(between(
        .data[[observed]],
        .data[[quant_low[11]]],
        .data[[quant_high[11]]]
      )),
      coverage60 = mean(between(
        .data[[observed]],
        .data[[quant_low[12]]],
        .data[[quant_high[12]]]
      )),
      coverage65 = mean(between(
        .data[[observed]],
        .data[[quant_low[13]]],
        .data[[quant_high[13]]]
      )),
      coverage70 = mean(between(
        .data[[observed]],
        .data[[quant_low[14]]],
        .data[[quant_high[14]]]
      )),
      coverage75 = mean(between(
        .data[[observed]],
        .data[[quant_low[15]]],
        .data[[quant_high[15]]]
      )),
      coverage80 = mean(between(
        .data[[observed]],
        .data[[quant_low[16]]],
        .data[[quant_high[16]]]
      )),
      coverage85 = mean(between(
        .data[[observed]],
        .data[[quant_low[17]]],
        .data[[quant_high[17]]]
      )),
      coverage90 = mean(between(
        .data[[observed]],
        .data[[quant_low[18]]],
        .data[[quant_high[18]]]
      )),
      coverage95 = mean(between(
        .data[[observed]],
        .data[[quant_low[19]]],
        .data[[quant_high[19]]]
      )),
      coverage005 = mean(.data[[observed]] < `quantile_0.5%`)
    ) %>%
    mutate(
      model = model.name,
      response = response
    )

  # same calculation by hour
  hour.scores.model <- hour_stats %>%
    group_by(hour) %>%
    summarise(
      rmse = sqrt(mean((.data[[observed]] - mean)^2)),
      across(!!stats.vec, mean),
      # crps = mean(crps),
      coverage05 = mean(between(
        .data[[observed]],
        .data[[quant_low[1]]],
        .data[[quant_high[1]]]
      )),
      coverage10 = mean(between(
        .data[[observed]],
        .data[[quant_low[2]]],
        .data[[quant_high[2]]]
      )),
      coverage15 = mean(between(
        .data[[observed]],
        .data[[quant_low[3]]],
        .data[[quant_high[3]]]
      )),
      coverage20 = mean(between(
        .data[[observed]],
        .data[[quant_low[4]]],
        .data[[quant_high[4]]]
      )),
      coverage25 = mean(between(
        .data[[observed]],
        .data[[quant_low[5]]],
        .data[[quant_high[5]]]
      )),
      coverage30 = mean(between(
        .data[[observed]],
        .data[[quant_low[6]]],
        .data[[quant_high[6]]]
      )),
      coverage35 = mean(between(
        .data[[observed]],
        .data[[quant_low[7]]],
        .data[[quant_high[7]]]
      )),
      coverage40 = mean(between(
        .data[[observed]],
        .data[[quant_low[8]]],
        .data[[quant_high[8]]]
      )),
      coverage45 = mean(between(
        .data[[observed]],
        .data[[quant_low[9]]],
        .data[[quant_high[9]]]
      )),
      coverage50 = mean(between(
        .data[[observed]],
        .data[[quant_low[10]]],
        .data[[quant_high[10]]]
      )),
      coverage55 = mean(between(
        .data[[observed]],
        .data[[quant_low[11]]],
        .data[[quant_high[11]]]
      )),
      coverage60 = mean(between(
        .data[[observed]],
        .data[[quant_low[12]]],
        .data[[quant_high[12]]]
      )),
      coverage65 = mean(between(
        .data[[observed]],
        .data[[quant_low[13]]],
        .data[[quant_high[13]]]
      )),
      coverage70 = mean(between(
        .data[[observed]],
        .data[[quant_low[14]]],
        .data[[quant_high[14]]]
      )),
      coverage75 = mean(between(
        .data[[observed]],
        .data[[quant_low[15]]],
        .data[[quant_high[15]]]
      )),
      coverage80 = mean(between(
        .data[[observed]],
        .data[[quant_low[16]]],
        .data[[quant_high[16]]]
      )),
      coverage85 = mean(between(
        .data[[observed]],
        .data[[quant_low[17]]],
        .data[[quant_high[17]]]
      )),
      coverage90 = mean(between(
        .data[[observed]],
        .data[[quant_low[18]]],
        .data[[quant_high[18]]]
      )),
      coverage95 = mean(between(
        .data[[observed]],
        .data[[quant_low[19]]],
        .data[[quant_high[19]]]
      )),
      coverage005 = mean(.data[[observed]] < `quantile_0.5%`)
    ) %>%
    mutate(
      model = model.name,
      response = response
    )

  # get global mean scores
  mean.scores <- daily.scores.model %>%
    # do.call(bind_rows,.) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
    mutate(
      model = model.name,
      response = response
    )

  return(list(
    global = mean.scores,
    day = daily.scores.model,
    # hour_stats = hour_stats,
    hour = hour.scores.model
  ))
}

save_samples <- function(
  t1,
  n.days = 213,
  quantile.seq = c(0.005, 0.025, 0.5, 0.975),
  n.samples = 1000,
  data,
  response,
  inla.object,
  path.fig,
  show.fig = FALSE,
  save.fig = TRUE,
  path.samples,
  mod.code,
  seq_filter = NULL
) {
  # day sequence to predict day-ahead power
  time_seq <- seq(
    from = t1,
    to = t1 + n.days * 24 * 3600,
    by = "day"
  )

  if (!is.null(seq_filter)) {
    time_seq <- time_seq[seq_filter]
  }

  # record initial time
  start_time <- Sys.time()

  # get model where each day in time_seq is predicted
  result <- lapply(
    time_seq,
    \(time0) {
      # re fit model
      new_fit <- fit_a_date(
        timet = time0,
        h = 24,
        dat.power = data,
        response = response,
        inla.object = inla.object,
        # verbose = TRUE
      )

      # labels according to model response
      resp.lab <- case_when(
        response == "actuals.cf" ~ "Wind Generation",
        response == "err.cf" ~ "Forecast error"
      )

      # save samples plot and simulate
      sim.obj <- simulation.plots.inla2(
        inla.model = new_fit,
        data = data,
        response = inla.object$.args$formula[[2]] %>% as.character(),
        t1 = time0,
        quSeq = quantile.seq,
        family = inla.object$.args$family,
        resp.lab = resp.lab,
        # ylim = c(0,1.02),
        nsamp = n.samples,
        show.fig = show.fig,
        save.fig = save.fig,
        path = path.fig,
        run.name = paste0(mod.code, time0)
        # sample.df = sample.test$samples,
        # legend.position = "bottom"
      )

      # path to store samples
      sample.path <- path.samples

      # create directory in case
      if (!dir.exists(sample.path)) {
        dir.create(sample.path, recursive = TRUE)
      }

      # write samples in a csv
      sim.obj$samples %>%
        as.data.frame() %>%
        setNames(paste0(
          "t_",
          seq(time0, time0 + 23 * 3600, by = "hour") %>%
            format(., "%Y-%m-%d_%H")
        )) %>%
        write.csv(
          .,
          file.path(sample.path, paste0(mod.code, time0, ".csv")),
          row.names = FALSE
        )
    }
  )

  end_time <- Sys.time()
  run_time <- end_time - start_time
  print(run_time)
}


ramps_model <- function(
  sample.path = file.path("~/Documents/proj2", "sample"),
  model.name,
  t1,
  ramp.thresholds = c(1, 2, 5, 10, 20) / 100,
  file_ext = ".csv",
  clipping = FALSE
  # response
) {
  # get response from model code
  # response = if_else(
  #   str_detect(model.name, "r_act_f"),
  #   "actuals.cf",
  #   if_else(str_detect(model.name, "r_err_f"), "err.cf", NA_character_)
  # )
  response <- case_when(
    str_detect(model.name, "actuals.cf|actuals-cf") ~ "actuals.cf",
    str_detect(model.name, "actuals") ~ "actuals",
    str_detect(model.name, "err.cf|err-cf") ~ "err.cf",
    str_detect(model.name, "err") ~ "err",
    TRUE ~ NA
  )
  #
  # capacity <- data.scaled[which(data.scaled$time == t1),] %>% pull("capacity")
  # forecast <- data.scaled[which(data.scaled$time == t1),] %>% pull("forecast.cf_orig")
  #
  # capacaity_trans <- case_when(
  #   str_detect(response,"actuals.cf")~ c(1,0),
  #   str_detect(response,"actuals")~ c(capacity,0),
  #   str_detect(response,"err.cf") ~ c(1,forecast),
  #   str_detect(response,"err") ~ c(capacity,forecast),
  #   TRUE ~ c(1,0)
  # )
  ramp.thresholds = c(1, 2, 5, 10, 20) / 100
  raw_samples <- file.path(
    sample.path,
    paste0(model.name, date(t1), file_ext)
  ) %>%
    fread(.)

  samp.tbl <- sample_rescale(
    t1,
    data.scaled,
    raw_samples,
    clipping,
    response
  )
  # Function to calculate percentage below thresholds
  percent_below_thresholds <- function(diff_values, thresholds) {
    sapply(thresholds, function(threshold) {
      mean(diff_values < threshold, na.rm = TRUE) * 100
    })
  }
  # browser()
  # read simulations
  scen.tbl <- #file.path(
    # sample.path,
    # paste0(model.name,date(t1),file_ext)) %>%
    # fread(.) %>%
    samp.tbl %>%
    mutate(id = 1:n()) %>%
    # select(-`NA`) %>%
    pivot_longer(cols = -id, names_to = "time", values_to = "value") %>%
    mutate(
      time = if_else(
        str_detect(time, "^NA$"), # Check if the value is "NA"
        NA_POSIXct_, # Keep it as NA if it's "NA"
        as.POSIXct(str_remove(time, "^t_"), format = "%Y-%m-%d_%H", tz = "PST")
      ),
      date = date(time)
    ) %>%
    group_by(id, date) %>%
    mutate(diff = abs(c(NA, diff(value)))) %>%
    reframe(
      across(
        .cols = diff,
        .fns = list(
          perc_below = ~ percent_below_thresholds(., ramp.thresholds)
        )
      )
    ) %>%
    mutate(threshold = rep(ramp.thresholds, n() / length(ramp.thresholds)))
  # pivot_wider(names_from = "id", values_from = response) %>%
  # as.data.frame()

  # compare with actuals - forecast error
  # fcst_dates <- seq(
  #   from = t1 + 0 * 60 * 60,
  #   to = t1 + 23 * 60 * 60, by = "hour")

  return(scen.tbl)
}

sample_rescale <- function(
  t1,
  data,
  samples,
  clipping = FALSE,
  response
) {
  # filtering data
  fcst_dates <- seq(
    from = t1 + 0 * 60 * 60,
    to = t1 + 23 * 60 * 60,
    by = "hour"
  )

  subsample <- data %>%
    filter(time %in% fcst_dates)

  # browser()
  additive_factor <- case_when(
    grepl("err.cf$", response) ~ subsample[["forecast.cf_orig"]],
    grepl("err$", response) ~ subsample[["forecast"]],
    TRUE ~ 0,
  )

  mult_factor <- case_when(
    grepl("actuals$", response) ~ subsample[["capacity"]],
    grepl("err$", response) ~ subsample[["capacity"]],
    TRUE ~ 1,
  )
  # rescaling
  # browser()
  rescaled <- apply(
    samples[, 1:24],
    1,
    \(row) {
      rescaled_sample <- (row + additive_factor) / mult_factor
      if (clipping) {
        rescaled_sample <- pmin(1, pmax(0, rescaled_sample))
      }
      rescaled_sample
    }
  ) %>%
    t() %>%
    as.data.frame() %>%
    setNames(format(fcst_dates, format = "t_%Y-%m-%d_%H"))

  # %>%
  # mutate(time = fcst_dates)

  rescaled
}

plot_nh_ahead <- function(
  data, # Input dataset
  mycols = c("actuals.cf", "forecast.cf"), # series to plot
  col_labels = c("observed", "forecast"),
  t0, # Starting time (t0)
  h = 24,
  show.fig = TRUE,
  legend.position = "bottom",
  position = c(0.2, 0.7),
  ylim = c(0, 1),
  xvar = "time",
  axis_2 = FALSE
) {
  # Number of hours ahead to plot (default = 24)

  # Subset the data based on selected columns and time window
  data.subset <- data %>%
    select(any_of("time"), any_of(c(xvar, mycols))) %>% # Ensure 'time' is always selected for plotting
    filter(.data[[xvar]] >= t0, .data[[xvar]] <= t0 + hours(h)) #%>%  # Filter time range from t0 to t0 + h hours
  # mutate(hour.lead = 1:n())         # Create an index for each time step (hourly lead)

  if (axis_2) {
    rng1 <- range(data.subset[[mycols[1]]], na.rm = TRUE)
    rng2 <- range(data.subset[[mycols[2]]], na.rm = TRUE)

    a <- diff(rng2) / diff(rng1)
    b <- rng2[1] - a * rng1[1]

    data.subset$y2_scaled <- (data.subset[[mycols[2]]] - b) / a
  } else {
    data.subset$y2_scaled <- data.subset[[mycols[2]]]
  }

  p <- ggplot(data.subset, aes(x = .data[[xvar]])) +
    geom_line(
      aes(
        y = y2_scaled,
        col = col_labels[2]
      )
    ) +
    geom_line(
      aes(
        y = !!sym(mycols[1]),
        col = col_labels[1]
      )
    ) +
    labs(
      col = "",
      y = "Normalised wind energy"
    ) +
    theme_bw() +
    theme(
      legend.position = legend.position,
      legend.position.inside = position,
      legend.background = element_rect(fill = NA, color = NA),
      legend.key = element_rect(fill = NA, color = NA)
    ) +
    coord_cartesian(ylim = ylim) +
    scale_color_lancet()

  if (axis_2) {
    p <- p +
      scale_y_continuous(
        sec.axis = sec_axis(
          ~ . * a + b,
          name = col_labels[2]
        )
      )
  }

  # # Initialize the ggplot object with the selected 'time' variable
  # fig <- ggplot(data.subset, aes(x = .data[[xvar]])) +

  #   # plot 2nd variable
  #   geom_line(aes(y = !!sym(mycols[2]), col = col_labels[2])) +
  #   # plot first variable as actuals
  #   geom_line(aes(y = !!sym(mycols[1]), col = col_labels[1])) +

  #   # Set plot limits, labels, theme, and color scale
  #   # coord_cartesian(ylim = c(0, 1)) +   # Limit the y-axis to [0,1]
  #   labs(col = "", y = "Normalised wind energy") + # Remove legend title
  #   theme_bw() + # Set a clean theme
  #   theme(
  #     legend.position = legend.position,
  #     legend.position.inside = position,
  #     legend.background = element_rect(fill = NA, color = NA),
  #     legend.key = element_rect(fill = NA, color = NA)
  #   ) + # Position legend at the inside
  #   coord_cartesian(ylim = ylim) +
  #   scale_color_lancet() # Apply Lancet color scheme

  # Print the plot
  if (show.fig) {
    print(p)
  }

  # Return the filtered data invisibly
  invisible(list(data = data.subset, plot = p))
}

compare_inla_models <- function(models) {
  # Ensure the input is a named list
  if (is.null(names(models))) {
    stop("The input must be a named list of models.")
  }

  # Extract WAIC and Mean Log-CPO for each model
  results <- lapply(names(models), function(model_name) {
    model <- models[[model_name]]

    waic <- if (!is.null(model$waic)) model$waic$waic else NA
    mean_log_cpo <- if (!is.null(model$cpo)) {
      mean(log(model$cpo$cpo), na.rm = TRUE)
    } else {
      NA
    }

    data.frame(
      Model = model_name,
      WAIC = waic,
      Mean_Log_CPO = mean_log_cpo
    )
  })

  # Combine into a single data frame
  results_df <- do.call(rbind, results)
}

extract_ar_details <- function(formula_string) {
  # Regex patterns
  ar1_pattern <- "f\\([^,]+,\\s*model\\s*=\\s*\"ar1\""
  ar_pattern <- "f\\([^,]+,\\s*model\\s*=\\s*\"ar\",.*?order\\s*=\\s*(\\d+)"
  # browser()
  # Check for ar1
  if (
    grepl(ar1_pattern, formula_string %>% paste(., collapse = " "), perl = TRUE)
  ) {
    return(1) # ar1 defaults to order = 1
  }

  # Check for ar and extract the order number
  ar_match <- regmatches(
    formula_string,
    regexpr(ar_pattern, formula_string, perl = TRUE)
  )
  if (length(ar_match) > 0 && nchar(ar_match) > 0) {
    order_match <- regmatches(ar_match, regexpr("\\d+", ar_match, perl = TRUE))
    return(as.integer(order_match))
  }

  # If neither ar1 nor ar is found, return 0
  return(0)
}

get_mod_stats <- function(
  fname,
  path = "~/Documents/proj2/model_objects",
  mod.obj = NULL
) {
  if (is.null(mod.obj)) {
    # read model
    mod.obj <- readRDS(file.path(path, fname))
  }

  # browser()

  features_vec <- c(
    "ws.w_group",
    "fcst_group",
    "fd_group",
    "eta deriv",
    "month",
    "hour",
    "ar1",
    "ar2",
    "t"
  )

  features_catalog <- data.frame(
    code = features_vec,
    readable = c(
      "wind speed",
      "forecast",
      "$\\Delta$ forecast",
      "\\Delta \\eta$",
      "Month",
      "Hour",
      "AR(1)",
      "AR(2)",
      "AR"
    )
  )
  readable_features <- data.frame(
    code = names(mod.obj$summary.random)
  ) %>%
    left_join(
      features_catalog,
      by = "code"
    )

  # check if eta deriv was included
  if (length(mod.obj$.args$family) == 1) {
    # extract info
    stats = data.frame(
      family = mod.obj$.args$family,
      response.code = mod.obj$.args$formula[[2]] %>% as.character(),
      fixed = paste(rownames(mod.obj$summary.fixed), collapse = "\n"),
      effects = paste(names(mod.obj$summary.random), collapse = "\n"),
      readable = paste(readable_features$readable, collapse = "\n"),
      m_rand = paste(mod.obj$model.random, collapse = "\n"),
      ar_order = mod.obj$.args$formula[[3]] %>% extract_ar_details(),
      waic = if (!is.null(mod.obj$waic)) mod.obj$waic$waic else NA,
      mean_log_cpo = if (!is.null(mod.obj$cpo)) {
        mean(log(mod.obj$cpo$cpo), na.rm = TRUE)
      } else {
        NA
      },
      mean_log_gcpo = if (!is.null(mod.obj$gcpo$gcpo)) {
        mean(log(mod.obj$gcpo$gcpo), na.rm = TRUE)
      } else {
        NA
      },
      etaderiv = FALSE
    )
  } else {
    # 2 likelihood models
    n = length(mod.obj$.args$data$eta) / 2

    stats = data.frame(
      family = mod.obj$.args$family[2],
      response.code = mod.obj$.args$formula[[2]] %>% as.character(),
      fixed = paste(rownames(mod.obj$summary.fixed), collapse = "\n"),
      effects = paste(names(mod.obj$summary.random), collapse = "\n"),
      m_rand = paste(mod.obj$model.random, collapse = "\n"),
      ar_order = mod.obj$.args$formula[[3]] %>% extract_ar_details(),
      waic = if (!is.null(mod.obj$waic)) {
        sum(mod.obj$waic$local.waic[-c(1:n)], na.rm = TRUE)
      } else {
        NA
      },
      mean_log_cpo = if (!is.null(mod.obj$cpo)) {
        mean(log(mod.obj$cpo$cpo[-c(1:n)]), na.rm = TRUE)
      } else {
        NA
      },
      mean_log_gcpo = if (!is.null(mod.obj$gcpo$gcpo)) {
        mean(log(mod.obj$gcpo$gcpo[-c(1:n)]), na.rm = TRUE)
      } else {
        NA
      },
      etaderiv = TRUE
    )
  }

  return(stats)
}


mask_inla <- function(
  data,
  columns = c("actuals.cf", "actuals", "err.cf", "err"),
  tt,
  h = 24 # hours
) {
  result <- data %>%
    mutate(across(any_of(columns), \(x) ifelse(time <= tt, x, NA)))
  return(result)
}

history_window <- function(
  scen.data,
  t,
  h = 24,
  window = 3, # months
  units = "months",
  mask = TRUE,
  ...
) {
  # browser()
  t0 <- t - months(window)
  t1 <- t + hours(h)

  # resulta <- scen.data %>%
  #   filter(time >= t0)

  # resultb <- resulta %>%
  #   filter(time <= t1)

  result <- scen.data %>% #resultb #%>%
    filter(time >= t0, time <= t1) %>%
    {
      if (mask) mask_inla(., tt = t) else .
    }
  # mutate(across(c(actuals.cf,actuals,err.cf), \(x) ifelse(time < t, x , NA)))
  return(result)
}


logit <- function(x) log(x / (1 - x))
expit <- function(x) exp(x) / (1 + exp(x))
trunc.cf <- function(x) pmax(pmin(1, x), 0)


build_effect <- function(features_opt) {
  # browser()
  # Convert options to a named list
  values <- as.list(features_opt)
  # Filter out NA values
  values <- values[!sapply(values, is.na)]
  # Construct key-value pairs for remaining arguments
  args <- paste(names(values), "=", unlist(values), collapse = ", ")
  # Create the final string
  result <- sprintf("f(%s)", args) %>%
    str_replace(., "feature =", "")
  return(result)
}


# rhs_formula <- function(features_vec){
#
#   # apply effect rules
#   features_opt <- data.frame(
#     feature0 = features_vec
#   ) %>%
#     mutate(
#       feature = case_when(
#         grepl("^ar",feature0) ~ "t",
#         TRUE ~ feature0
#       ),
#       model = case_when(
#         grepl("^ar",feature0) ~ ifelse(substr(feature0,3,3)=="1","'ar1'","'ar'"),
#         TRUE ~ "'rw2'"
#       ),
#       order = case_when(
#         model == "'ar'" ~ 2,
#         TRUE ~ NA
#       ),
#       cyclic = case_when(
#         feature %in% c("month","hour") ~ TRUE,
#         TRUE ~ NA
#       ),
#       hyper = case_when(
#         model == "'rw2'" ~ "hyper.rw2",
#         model == "'ar1'" ~ "hyper.ar1",
#         model == "'ar'" ~ "hyper.ar2",
#         TRUE ~ NA
#       )
#     )
#
#   # convert to effects syntax
#   effects <- apply(features_opt %>% select(-feature0), 1, build_effect)
#
#   # terms for rhs
#   rhs <- paste0(effects, collapse = " + ")
#   # browser()
#   return(rhs)
# }

rhs_formula <- function(features_vec) {
  # browser()
  # check if eta derivative will be used
  if ("etaderiv" %in% features_vec) {
    # set status to add eta model at the end
    eta = TRUE
    # remove etaderiv from features vec
    features_vec <- features_vec[-which(features_vec == "etaderiv")]
  } else {
    eta = FALSE
  }

  # apply effect rules
  features_opt <- data.frame(
    feature0 = features_vec # features code names
  ) %>%
    # construction effects for each type of feature
    mutate(
      # relabelling AR terms to time index
      feature = case_when(
        grepl("^ar", feature0) ~ "t",
        TRUE ~ feature0
      ),
      # model for each feature
      model = case_when(
        grepl("^ar", feature0) ~ case_when(
          substr(feature0, 3, 3) == "1" ~ "'ar1'",
          substr(feature0, 3, 3) == "2" ~ "'ar'"
        ),
        TRUE ~ "'rw2'"
      ),
      # order for AR terms
      order = case_when(
        model == "'ar'" ~ 2,
        TRUE ~ NA
      ),
      # cyclic option for time effects
      cyclic = case_when(
        feature %in% c("month", "hour") ~ TRUE,
        TRUE ~ NA
      ),
      # priors for each model
      hyper = case_when(
        grepl("ar1g|ar2g", feature0) ~ "NULL",
        model == "'rw2'" ~ "hyper.rw2",
        model == "'ar1'" ~ "hyper.ar1",
        model == "'ar'" ~ "hyper.ar2",
        TRUE ~ NA
      ),
      group = case_when(
        grepl("ar1g$|ar2g$", feature0) ~ "site_id",
        TRUE ~ NA
      ),
      values = case_when(
        grepl("ar1g$|ar2g$", feature0) ~ "sort(unique(t))",
        TRUE ~ NA
      )
    )

  # convert to effects syntax
  effects_vec <- apply(features_opt %>% select(-feature0), 1, build_effect)

  # finish with eta related effects if required
  if (eta) {
    eta_effects <- c(
      "f(eta, w, model = 'iid',
        hyper = list(prec = list(initial = -15, fixed = TRUE)))", # linear predictor with fixed precision
      "f(eta.1, copy = 'eta',
        hyper = list(beta = list(initial = -0.1, fixed = FALSE)),
        range = c(-0.4,0.4))", # effect for positive part of deriv
      "f(eta.2, w2, copy = 'eta',  same.as = 'eta.1')" # effect for negative part of deriv
    )
    effects_vec <- c(effects_vec, eta_effects)
  }

  # terms for rhs
  rhs <- paste0(effects_vec, collapse = " + \n")
  # browser()
  return(rhs)
}

# build_formula <- function(
#     model_type,
#     features_vec
# ){
#
#   formula <- as.formula(
#     paste0(
#       model_type$response,
#       " ~ ",
#       rhs_formula(features_vec)
#     )
#   )
#   return(formula)
# }

build_formula <- function(
  model_type,
  features_vec
) {
  # building left hand side
  lhs <- model_type$response
  # renaming response for multiple likelihood models
  if ("etaderiv" %in% features_vec) {
    lhs <- paste0("Y_", lhs)
  }

  # build formula
  formula <- as.formula(
    paste0(
      lhs, # LHS
      " ~ ",
      rhs_formula(features_vec) # RHS
    )
  )
  return(formula)
}

# fit_inla_model <- function(
#     model_type,
#     features_vec,
#     data,
#     ini.theta,
#     restart = TRUE,
#     cens = 0.001,
#     verbose = FALSE,
#     mycontrol.inla = list(),
#     gcpo = TRUE,
#     n.groups = 3,
#     ...
# ){
#   # family_opts <- ifelse(
#   #   model_type$family == "beta",
#   #   list(beta.censor.value = cens,
#   #        control.link = list(model = "default")),
#   #   NULL
#   # )
#   family_opts <- list(beta.censor.value = cens,
#                       control.link = list(model = "default"))
#   # browser()
#   # exclude zeroes
#   if(model_type$family %in% c("weibull","gamma")){
#     data <- data %>%
#       filter(.data[[model_type$response]]>0)
#   }
#
#   # build initials from arguments
#   mode_opts <- if (identical(ini.theta, list())) {
#     ini.theta
#   } else {
#     list(theta = ini.theta, restart = restart)
#   }
#
#   # override gcpo option for beta models
#   if(model_type$family %in% c("beta")){
#     gcpo = FALSE
#   }
#
#   # build gcpo options
#   gcpo_opts <- if (gcpo){
#     control.gcpo(enable = TRUE, num.level.sets = n.groups)
#   } else {
#     control.gcpo()
#   }
#
#   # if(model_type$family == "beta"){
#   #   family_opts
#   # }
#   # browser()
#   model_formula <- build_formula(model_type, features_vec)
#   cat("Fitting formula:\n")
#   print(model_formula)
#   cat(sprintf("Using %s family", model_type$family))
#
#   return(
#     inla(
#       formula = model_formula,
#       data = data,
#       family = model_type$family,
#       control.mode = mode_opts,
#       control.family = family_opts,
#       control.compute = list(
#         config = TRUE, waic = TRUE, cpo = TRUE,
#         control.gcpo = gcpo_opts),
#       control.predictor = list(compute = TRUE, link = 1),
#       control.inla = mycontrol.inla,
#       verbose = verbose
#     )
#   )
# }

fit_inla_model <- function(
  model_type,
  features_vec,
  data0,
  ini.theta,
  restart = TRUE,
  cens = 0.001,
  verbose = FALSE,
  mycontrol.inla = list(),
  gcpo = TRUE,
  n.groups = 3,
  ...
) {
  # family_opts <- ifelse(
  #   model_type$family == "beta",
  #   list(beta.censor.value = cens,
  #        control.link = list(model = "default")),
  #   NULL
  # )

  # browser()
  # exclude zeroes
  if (model_type$family %in% c("weibull", "gamma")) {
    data0 <- data0 %>%
      filter(
        .data[[model_type$response]] > 0 | is.na(.data[[model_type$response]])
      )
  }

  # restructure data for eta deriv type
  if ("etaderiv" %in% features_vec) {
    n <- nrow(data0)
    # change response
    updated_response <- paste0("Y_", model_type$response)

    # change data
    data <- with(
      data0,
      list(
        # Dynamically set the name using updated_response
        intercept = c(rep(1, n), rep(NA, n)), # 1st likelihood intercept, nothing for 2nd as it's stored in eta
        # wind = c(ws.w, rep(NA, n)), # x data goes in 1st likelihood
        ws.w_group = c(ws.w_group, rep(NA, n)),
        fcst_group = c(fcst_group, rep(NA, n)),
        fd_group = c(fd_group, rep(NA, n)),
        t = c(t, rep(NA, n)),
        time = c(time, rep(NA, n)),
        month = c(month, rep(NA, n)),
        hour = c(hour, rep(NA, n)),
        site_id = c(site_id, rep(NA, n)),
        eta = c(1:n, 1:n), # eta indices
        w = c(rep(-1, n), rep(1, n)), # weights for eta effect: -1 to copy lin.predictor, 1 to be part of likelihood
        eta.1 = c(rep(NA, n), 1:n), # indices for positive part of derivative
        eta.2 = c(rep(NA, n), NA, 1:(n - 1)), # shifted indices for negative part of derivative
        w2 = c(rep(0, n), rep(-1, n))
      ) # weights for negative part of derivative
    )
    # adding properly named response
    data[[updated_response]] <- cbind(
      c(rep(0, n), rep(NA, n)), # fake zeros (1st likelihood)
      c(rep(NA, n), data0[[model_type$response]])
    )

    # change family
    family <- c("gaussian", model_type$family)

    # control family
    ctrl_family <- list(
      list(hyper = list(prec = list(initial = 15, fixed = TRUE))), # 1st lik fixed precision
      list(
        #hyper = list(theta = list(initial = log(1/s^2), fixed = FALSE)),
        control.link = list(model = "default"),
        beta.censor.value = cens
      ) # 2nd beta lik non-fixed precision
    )
  } else {
    data <- data0
    family <- model_type$family
    ctrl_family <- list(
      beta.censor.value = cens,
      control.link = list(model = "default")
    )
  }

  # build initials from arguments
  mode_opts <- if (identical(ini.theta, list())) {
    ini.theta
  } else {
    list(theta = ini.theta, restart = restart)
  }

  # override gcpo option for beta models
  if (model_type$family %in% c("beta")) {
    gcpo = FALSE
  }

  # build gcpo options
  gcpo_opts <- if (gcpo) {
    control.gcpo(enable = TRUE, num.level.sets = n.groups)
  } else {
    control.gcpo()
  }

  # if(model_type$family == "beta"){
  #   family_opts
  # }

  model_formula <- build_formula(model_type, features_vec)
  cat("Fitting formula:\n")
  print(model_formula)
  cat(sprintf("Using %s family\n", model_type$family))
  # browser()
  return(
    inla(
      formula = model_formula,
      data = data,
      family = family,
      control.mode = mode_opts,
      control.family = ctrl_family,
      control.compute = list(
        config = TRUE,
        waic = TRUE,
        cpo = TRUE,
        control.gcpo = gcpo_opts
      ),
      control.predictor = list(compute = TRUE, link = 1),
      control.inla = mycontrol.inla,
      verbose = verbose
    )
  )
}

DIC <- function(x) {
  data.frame(
    mean.deviance = x$dic$mean.deviance,
    p.eff = x$dic$p.eff,
    dic = x$dic$dic,
    waic = x$waic$waic
  )
}

CV <- function(x, response) {
  LS <- -mean(log(x$cv))

  n <- length(x$cv)
  # Computation of E[Y_i | y_{I_i}] using trapezoid numerical integration #
  expectation = numeric(n)
  for (i in 1:n) {
    mu = x$mean[i]
    sd = x$sd[i]
    xx = seq(mu - 6 * sd, mu + 6 * sd, length.out = 100)
    yy = xx * dnorm(xx, mean = mu, sd = sd)
    expectation[i] = pracma::trapz(xx, yy)
  }
  MSPE <- mean((expectation - response)^2)

  data.frame(LS = LS, MSPE = MSPE)
}

extract_score_model <- function(mod.obj) {
  # browser()
  # check if model has two likelihoods
  # n <- if (length(mod.obj$.args$family) > 1) length(mod.obj$.args$data$intercept) / 2 else NULL

  if (length(mod.obj$.args$family) > 1) {
    n <- length(mod.obj$.args$data$intercept) / 2
    result <- data.frame(
      waic = if (!is.null(mod.obj$waic)) {
        sum(mod.obj$waic$local.waic[-c(1:n)], na.rm = TRUE)
      } else {
        NA
      },
      mean_log_cpo = if (!is.null(mod.obj$cpo$cpo)) {
        mean(log(mod.obj$cpo$cpo[-c(1:n)]), na.rm = TRUE)
      } else {
        NA
      },
      mean_log_gcpo = if (!is.null(mod.obj$gcpo$gcpo)) {
        mean(log(mod.obj$gcpo$gcpo[-c(1:n)]), na.rm = TRUE)
      } else {
        NA
      }
    )
  } else {
    result <- data.frame(
      waic = if (!is.null(mod.obj$waic)) mod.obj$waic$waic else NA,
      mean_log_cpo = if (!is.null(mod.obj$cpo)) {
        mean(log(mod.obj$cpo$cpo), na.rm = TRUE)
      } else {
        NA
      },
      mean_log_gcpo = if (!is.null(mod.obj$gcpo$gcpo)) {
        mean(log(mod.obj$gcpo$gcpo), na.rm = TRUE)
      } else {
        NA
      }
    )
  }
  return(result)
}


get_regime_models <- function(
  model_file,
  in_path,
  data,
  group_var = "power.group",
  out_path = "~/Documents/proj2/model_objects/regime"
) {
  # browser()
  p_groups <- data[[group_var]] %>% unique()

  for (k in p_groups) {
    resp_columns = c("actuals.cf", "actuals", "err.cf", "err")

    # get data
    regime_data <- history_window(
      data,
      t1,
      window = 18,
      mask = T
    ) %>%
      # mask other regimes
      mutate(across(any_of(resp_columns), \(x) {
        ifelse(.data[[group_var]] %in% k, x, NA)
      }))

    # load previous model
    base_model <- readRDS(
      file.path(in_path, model_file)
    )

    base_model_stats <- get_mod_stats("", "", base_model)
    # browser()
    # attempt model train
    # regime_model <- #tryCatch(
    #   fit_a_date(
    #     timet = t1,
    #     h = 24,
    #     dat.power = regime_data,
    #     response = base_model_stats$response.code,
    #     inla.object = base_model,
    #     # verbose = TRUE
    #   )
    # )
    regime_model <- tryCatch(
      fit_a_date(
        timet = t1,
        h = 24,
        dat.power = regime_data,
        response = base_model_stats$response.code,
        inla.object = base_model,
        # verbose = TRUE
      ),
      error = function(e) {
        message("Error in model fitting: ", e$message)
        return(NULL)
      }
    )
    # save result if successful
    # browser()
    regmod_name <- sprintf(
      model_file %>%
        sub("(\\.rds)$", "_reg_%s\\1", .),
      substr(k, 1, 1)
    )

    if (!is.null(regime_model)) {
      saveRDS(
        regime_model,
        file.path(
          out_path,
          regmod_name
        )
      )
    }
  }
}


hist_wdens <- function(data, response, family = "gaussian") {
  library(ggplot2)
  # library(MASS)  # for fitdistr()

  # Extract and clean values
  values <- data[[response]]
  values <- values[!is.na(values)]

  # Remove non-positive values for certain families
  if (family %in% c("gamma", "weibull") && any(values <= 0)) {
    prop_neg <- mean(values <= 0)
    message(sprintf(
      "Found %.2f%% non-positive values in '%s'. These were removed before fitting the %s distribution.",
      prop_neg * 100,
      response,
      family
    ))
    values <- values[values > 0]
  }

  # Beta-specific check: must be in (0, 1)
  if (family == "beta") {
    if (any(values <= 0 | values >= 1)) {
      prop_out <- mean(values <= 0 | values >= 1)
      message(sprintf(
        "Beta distribution only supports values strictly in (0,1). %.2f%% of values are outside this range.",
        prop_out * 100
      ))
    }
  }

  # Base plot using cleaned values
  p <- ggplot(data.frame(values = values), aes(x = values)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 50,
      fill = "grey80",
      color = "black"
    )

  # Add density overlay
  if (family == "gaussian") {
    mu <- mean(values)
    sigma <- sd(values)
    p <- p +
      stat_function(
        fun = dnorm,
        args = list(mean = mu, sd = sigma),
        color = "darkred",
        lwd = 1
      )
  } else if (family == "gamma") {
    shape <- (mean(values)^2) / var(values)
    rate <- mean(values) / var(values)
    p <- p +
      stat_function(
        fun = dgamma,
        args = list(shape = shape, rate = rate),
        color = "darkblue",
        lwd = 1
      )
  } else if (family == "weibull") {
    fit <- suppressWarnings(MASS::fitdistr(values, "weibull"))
    shape <- fit$estimate["shape"]
    scale <- fit$estimate["scale"]
    p <- p +
      stat_function(
        fun = dweibull,
        args = list(shape = shape, scale = scale),
        color = "forestgreen",
        lwd = 1
      )
  } else if (family == "beta") {
    # Method of moments estimation
    m <- mean(values)
    v <- var(values)
    alpha <- m * ((m * (1 - m)) / v - 1)
    beta <- (1 - m) * ((m * (1 - m)) / v - 1)
    p <- p +
      stat_function(
        fun = dbeta,
        args = list(shape1 = alpha, shape2 = beta),
        color = "purple",
        lwd = 1
      )
  } else {
    stop("Unsupported family. Use 'gaussian', 'gamma', 'weibull', or 'beta'.")
  }

  # Final touches
  p <- p +
    labs(
      title = paste("Histogram with", family, "Density Overlay"),
      x = response,
      y = "Density"
    )
  print(p)
}

qq_wdens <- function(data, response, family = "gaussian") {
  library(ggplot2)
  # library(MASS)  # for fitdistr()

  # Extract and clean values
  values <- data[[response]]
  values <- values[!is.na(values)]

  # Handle zero or negative values
  if (family %in% c("gamma", "weibull") && any(values <= 0)) {
    prop_neg <- mean(values <= 0)
    message(sprintf(
      "Found %.2f%% non-positive values in '%s'. These were removed before fitting the %s distribution.",
      prop_neg * 100,
      response,
      family
    ))
    values <- values[values > 0]
  }

  # Handle Beta-specific domain check
  if (family == "beta") {
    if (any(values <= 0 | values >= 1)) {
      prop_out <- mean(values <= 0 | values >= 1)
      message(sprintf(
        "Found %.2f%% values outside (0,1) in '%s'. These were removed before fitting the beta distribution.",
        prop_out * 100,
        response
      ))
      values <- values[values > 0 & values < 1]
    }
  }

  # Sort and prepare quantiles
  n <- length(values)
  empirical <- sort(values)
  probs <- ppoints(n)

  # Compute theoretical quantiles
  if (family == "gaussian") {
    mu <- mean(values)
    sigma <- sd(values)
    theoretical <- qnorm(probs, mean = mu, sd = sigma)
    dist_color <- "darkred"
  } else if (family == "gamma") {
    shape <- (mean(values)^2) / var(values)
    rate <- mean(values) / var(values)
    theoretical <- qgamma(probs, shape = shape, rate = rate)
    dist_color <- "darkblue"
  } else if (family == "weibull") {
    fit <- suppressWarnings(MASS::fitdistr(values, "weibull"))
    shape <- fit$estimate["shape"]
    scale <- fit$estimate["scale"]
    theoretical <- qweibull(probs, shape = shape, scale = scale)
    dist_color <- "forestgreen"
  } else if (family == "beta") {
    m <- mean(values)
    v <- var(values)
    alpha <- m * ((m * (1 - m)) / v - 1)
    beta <- (1 - m) * ((m * (1 - m)) / v - 1)
    theoretical <- qbeta(probs, shape1 = alpha, shape2 = beta)
    dist_color <- "purple"
  } else {
    stop("Unsupported family. Use 'gaussian', 'gamma', 'weibull', or 'beta'.")
  }

  # Q-Q plot
  p <- ggplot(
    data.frame(empirical = empirical, theoretical = theoretical),
    aes(x = theoretical, y = empirical)
  ) +
    geom_point(color = dist_color, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(
      title = paste(response, "Q-Q Plot against", family, "Distribution"),
      x = "Theoretical Quantiles",
      y = "Empirical Quantiles"
    ) +
    theme_minimal()

  print(p)
}
plot_reliability <- function(
  global_scores,
  model_list,
  my_palette = scale_color_lancet()(9),
  show.fig = TRUE,
  pgroup = FALSE,
  code_subset = NULL
) {
  ordered.levels <- c("NPB", "NPN", "NPG", "PN", "PG", "EN", "NEN", "nonpar.")
  plot.data <- global_scores %>%
    bind_rows() %>%
    {
      if (pgroup) select(., c(1:29, 31)) else select(., 1:29)
    } %>%
    # bind_cols(
    #   model_list %>% select(response,family, transformation)
    # ) %>%
    left_join(
      model_cat,
      by = "model"
    ) %>%
    bind_rows(nonpar.tbl) %>%
    pivot_longer(
      cols = coverage05:coverage95,
      names_to = "level",
      values_to = "empirical"
    ) %>%
    mutate(nominal = (substr(level, 9, 10) %>% as.numeric()) / 100) %>%
    mutate(
      label = sprintf(
        "%s %s %s",
        ifelse(transformation == "none", "unnorm.", "norm."),
        ifelse(grepl("act", response), "power", "error"),
        family
      ),
      label2 = ifelse(
        family != "nonparametric",
        sprintf(
          "%s%s%s",
          ifelse(transformation == "none", "", "N"),
          ifelse(grepl("act", response), "P", "E"),
          case_when(
            family == "gaussian" ~ "N",
            TRUE ~ toupper(substr(family, 1, 1))
          )
        ),
        "nonpar."
      ) %>%
        factor(., levels = ordered.levels)
    ) %>%
    {
      if (!is.null(code_subset)) {
        filter(., as.character(label2) %in% code_subset)
      } else {
        .
      }
    } %>%
    {
      if (pgroup) {
        mutate(
          .,
          label2 = paste(
            label2,
            power.group %>%
              gsub("20%$", "", .)
          ) %>%
            gsub(" NA$", "", .) %>%
            factor(
              .,
              levels = c(
                "NPB low",
                "NPB mid",
                "NPB high",
                "NEN low",
                "NEN mid",
                "NEN high",
                "nonpar."
              )
            )
        )
      } else {
        .
      }
    }
  # browser()
  # my_palette <- pal_lancet()(9) %>% rev()
  # names(my_palette) <- plot.data$label %>% unique()
  p <- plot.data %>%
    ggplot() +
    geom_abline(aes(slope = 1, intercept = 0), col = "darkgray") +
    geom_line(
      aes(
        nominal,
        empirical,
        # col = paste0(
        #   "response: ",response,", family: ", family,", " ,
        #   main_term %>% substr(.,1,3) %>% toupper())),
        col = label2
      ),
      # col = paste0(
      # "response: ",response,", family: ", family,", " )),
      # main_term %>% substr(.,1,3) %>% toupper())),
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
      legend.position.inside = c(.8, .3),
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
    # scale_color_lancet()
    scale_color_manual(values = my_palette)
  if (show.fig) {
    print(p)
  }
  invisible(p)
}


power_ramp_freq <- function(
  samples_path,
  model_name,
  version = 1,
  ...
) {
  date_seq <- list.files(samples_path) %>% # files available in path
    sub(".*_t(.*)\\.csv.*", "\\1", .) # get dates
  # browser()
  model_code <- sub("(.*_t).*", "\\1", list.files(samples_path)[1])

  model.ramp.df <- lapply(
    date_seq,
    \(t) {
      ramps_model(
        sample.path = samples_path,
        model.name = model_code,
        t1 = t %>% as.POSIXct(., tz = "PST"),
        file_ext = ".csv.gz",
        ...
      ) %>%
        ungroup()
    }
  ) %>%
    bind_rows() %>%
    group_by(threshold) %>%
    summarise(
      # l_025 = quantile(ramp.freq,0.025, na.rm = T),
      # l_975 = quantile(ramp.freq,0.975, na.rm = T),
      ramp.freq = mean(diff_perc_below, na.rm = T),
      # model = model.list[1],
      .groups = "drop"
    ) %>%
    mutate(model = model_name, version = version)
}
