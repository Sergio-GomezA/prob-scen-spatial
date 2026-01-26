require(mgcv)

mysplinederiv <- function(
  y,
  x = 1:length(y),
  binsize = length(x),
  wholevec = TRUE,
  lead = 1
) {
  b <- gam(y ~ s(x), family = gaussian(link = "log"))
  fd <- derivatives(b, n = binsize)$.derivative
  if (!wholevec) {
    midpoint <- floor(binsize / 2)
    result <- fd[midpoint + 0:2]
  } else {
    result <- fd
  }
  result
}

mysplinederiv <- function(
  y,
  x = 1:length(y),
  binsize = length(x),
  wholevec = TRUE,
  lead = 1,
  k = -1, # Default: let mgcv decide knots
  bs = "tp", # Basis type: thin-plate spline by default
  method = "REML", # Penalization method
  na_action = c("omit", "interpolate") # How to handle NA values
) {
  # Match the `na_action` argument
  na_action <- match.arg(na_action)

  if (anyNA(y)) {
    if (na_action == "omit") {
      # Remove NA values along with corresponding x values
      na_idx <- !is.na(y)
      y <- y[na_idx]
      x <- x[na_idx]
    } else if (na_action == "interpolate") {
      # Interpolate NA values
      y <- zoo::na.approx(y, x = x, na.rm = FALSE) # Ensure consistent length
    }
  }

  # Check for sufficient data
  if (length(y) < 3) {
    stop("Insufficient data for spline fitting.")
  }

  # Fit GAM with specified basis and penalization
  b <- gam(
    y ~ s(x, k = k, bs = bs),
    family = gaussian(link = "log"),
    method = method
  )

  # Compute derivatives
  fd <- derivatives(b, n = binsize)$.derivative

  if (!wholevec) {
    # Extract specific values if only a subset is needed
    midpoint <- floor(binsize / 2)
    result <- fd[midpoint + 0:2]
  } else {
    result <- fd
  }

  list(derivative = result, model = b)
}

dist.estimate <- function(x, y, s, lowr = -400) {
  pdf.est <- approxfun(x, y)
  cdf <- integrate(
    f = pdf.est,
    lower = lowr,
    upper = s,
    # x = x, y=y,
    stop.on.error = FALSE
  )$value

  # result <- list(pdf=pdf,cdf=cdf)
  return(cdf)
}

getepidens <- function(
  pat0,
  data,
  response = "err.cf",
  runquiet = T,
  remove.tails = T,
  probs = c(.01, .99),
  ...
) {
  # initilising epi estimation
  soft <- setup.softinfo(10, ...)

  # filter data
  dat.filtered <- data %>%
    filter(pattern %in% pat0) %>%
    {
      if (remove.tails) {
        filter(
          .,
          between(
            data[[response]],
            quantile(data[[response]], probs[1]),
            quantile(data[[response]], probs[2])
          )
        )
      } else {
        .
      }
    }

  if (runquiet) {
    epidensity <- quiet(expepi(
      data = dat.filtered[[response]], # only works on this column
      softinfo = soft,
      postproc.controls = postproc.control(pic.types = NULL)
    ))
  } else {
    epidensity <- expepi(
      data = dat.filtered[[response]], # only works on this column
      softinfo = soft,
      postproc.controls = postproc.control(pic.types = NULL)
    )
  }

  epidensity
}

# PENDING check if approx works better
epi.cdf <- function(s, epidens0) {
  pdf.est <- approxfun(epidens0$x.pts, epidens0$y.est)
  cdf <- integrate(
    f = pdf.est,
    lower = epidens0$epiparameters$m0,
    upper = s,
    # x = x, y=y,
    stop.on.error = FALSE
  )$value
  cdf <- min(max(cdf, 0), 1)
  # result <- list(pdf=pdf,cdf=cdf)
  return(cdf)
}


# Hypo-Distance between Patterns
hypo_distance <- function(distrib_f, distrib_g, x) {
  # initialize d=0
  d <- 0
  K <- length(x)
  for (k in 1:K) {
    if (distrib_f[k] > distrib_g[k]) {
      j0 <- which.min((x[k] - x)^2 + (distrib_f[k] - distrib_g)^2)
      mind_temp <- (x[k] - x[j0])^2 + (distrib_f[k] - distrib_g[j0])^2
      d <- max(d, mind_temp)
    } else {
      j0 <- which.min((x[k] - x)^2 + (distrib_g[k] - distrib_f)^2)
      mind_temp <- (x[k] - x[j0])^2 + (distrib_g[k] - distrib_f[j0])^2
      d <- max(d, mind_temp)
    }
  }
  return(d)
}


consistentepi <- function(epidens_A, epidens_B, numevalpts = 10000) {
  min_x <- min(epidens_A$epiparameters$m0, epidens_B$epiparameters$m0)
  max_x <- max(epidens_A$epiparameters$mN, epidens_B$epiparameters$mN)
  x <- seq(min_x, max_x, length.out = numevalpts)

  epidens_A_common <- sapply(
    x,
    \(x) approx(epidens_A$x.pts, epidens_A$y.est, x)$y
  )

  epidens_B_common <- sapply(
    x,
    \(x) approx(epidens_B$x.pts, epidens_B$y.est, x)$y
  )

  result <- list(densA = epidens_B_common, densB = epidens_B_common, x = x)
  result
}


consistentepi <- function(epidens_A, epidens_B, numevalpts = 10000) {
  min_x <- min(epidens_A$epiparameters$m0, epidens_B$epiparameters$m0)
  max_x <- max(epidens_A$epiparameters$mN, epidens_B$epiparameters$mN)
  x <- seq(min_x, max_x, length.out = numevalpts)

  epidens_A_common <- approx(epidens_A$x.pts, epidens_A$y.est, x)$y
  epidens_B_common <- approx(epidens_B$x.pts, epidens_B$y.est, x)$y

  result <- list(
    densA = replace(epidens_A_common, is.na(epidens_A_common), 0),
    densB = replace(epidens_B_common, is.na(epidens_B_common), 0),
    x = x
  )
  result
}

distance_epi <- function(epidensA, epidensB) {
  # get estimates on same values of x
  comparable <- consistentepi(epidensA, epidensB)

  d <- hypo_distance(comparable$densA, comparable$densB, comparable$x)
  d
}

clust_weights <- function(hypo_dist, quantile_cut = 0.25) {
  # weight_mat <- ifelse(
  #   hypodist > tau,
  #   0,
  #   -hypo_dist/max(hypo_dist)+1
  # )

  above_threshold <- as.matrix(hypo_dist) >
    quantile(as.matrix(hypo_dist), quantile_cut)

  weight_mat <- as.matrix(-hypo_dist / max(hypo_dist) + 1)

  weight_mat[above_threshold] <- 0
  weight_mat
}

regimelabel <- function(x, alpha = .40, data) {
  # empirical cdf
  empirical_x <- mean(x <= data, na.rm = T)

  if (empirical_x <= alpha | empirical_x > 1 - alpha) {
    # obtaining cut-points for low and high regime
    data_endpoints <- quantile(data, probs = c(alpha, 1 - alpha), na.rm = T)
  } else {
    # mid regime cut-points
    data_endpoints <- quantile(
      data,
      probs = c(empirical_x - alpha / 2, empirical_x + alpha / 2),
      na.rm = T
    )
  }

  if (x < data_endpoints[1]) {
    # if x is in the first alpha % of the data
    index <- which(data <= data_endpoints[1])
  } else {
    if (x > data_endpoints[2]) {
      # if x is in the last alpha % of the data
      index <- which(data > data_endpoints[2])
    } else {
      index <- which(
        data > data_endpoints[1] &
          data <= data_endpoints[2]
      )
    }
  }

  index
}


add_patterns <- function(
  data,
  sel.column = "forecast",
  method = c("smooth_deriv", "differences"),
  binsize = 48,
  probs = c(0, .3, .7, 1)
) {
  # Ensure method is a valid choice and set default to "smooth_deriv"
  method <- match.arg(method)

  N <- nrow(data)

  if (method == "smooth_deriv") {
    # obtain derivatives

    fcstderiv <- sapply(1:ceiling(N / binsize), \(x) {
      mysplinederiv(data[[sel.column]][
        ((x - 1) * binsize + 1):min(x * binsize, N)
      ])
    })

    # EENDING: find out which one is getting the warning

    # combine in single row
    fcstderiv <- do.call(c, fcstderiv)
  } else {
    if (method == "differences") {
      fcstderiv <- data %>%
        mutate(
          zoo::na.approx(.data[[sel.column]], na.rm = FALSE),
          change = c(NA, diff(.data[[sel.column]]))
        ) %>%
        pull(change)
    }
  }

  # get derivatives quantiles to partition in patterns
  deriv_cuts <- quantile(fcstderiv, probs, na.rm = TRUE)
  # plot(density(fcstderiv[between(fcstderiv,quantile(fcstderiv,.1),quantile(fcstderiv,.9))]),
  #      main = "Forecast Derivatives density (approximation through splines)")

  # label patterns for each point
  deriv_cat <- cut(
    fcstderiv,
    breaks = deriv_cuts,
    labels = c("-", "0", "+"),
    include.lowest = T
  )

  # join local pattern with pattern before and after
  pattern <- paste(
    deriv_cat[-c(N - 1, N)],
    deriv_cat[-c(1, N)],
    deriv_cat[-c(1, 2)]
  )

  # add to data set
  data$fcstderiv <- fcstderiv
  data$pattern <- (c(NA, pattern, NA))
  data
}


dist.custom.mc <- function(list, dist.fun, ...) {
  # get pairs
  my_pair <- combn(1:length(list), 2, simplify = FALSE)

  # calculate distance over pairs
  out <- mclapply(
    my_pair,
    \(pair) {
      distance_epi(epidensA = list[[pair[1]]], epidensB = list[[pair[2]]])
    },
    ...
  )

  # convert to distance matris
  hypodist <- matrix(0, length(list), length(list))
  hypodist[lower.tri(hypodist)] <- out %>% unlist()
  # add names
  colnames(hypodist) <- names(list)
  rownames(hypodist) <- names(list)
  hypodist <- as.dist(hypodist + t(hypodist), diag = T, upper = T)
  hypodist
}


get_clusters <- function(
  data,
  response = "err.cf",
  minimum.size = 20,
  mc = detectCores() - 2,
  ...
) {
  # quantify amount of data in each pattern
  pattern.size <- data %>%
    na.omit() %>%
    # mutate(pattern = factor(pattern)) %>%
    group_by(pattern) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(n) %>%
    mutate(size = ifelse(n < minimum.size, "small", "large"))

  # cluster patterns
  largep.list <- pattern.size$pattern[pattern.size$size == "large"]
  names(largep.list) <- largep.list

  # scale data for stability

  scaled.err <- data[[response]] %>% scale()
  err.ctr <- attr(scaled.err, which = "scaled:center")
  err.scl <- attr(scaled.err, which = "scaled:scale")

  data_fixed <- data %>%
    mutate(!!response := na.approx(.data[[response]])) %>% # interpolate NAs
    mutate(
      !!response := scale(
        .data[[response]],
        center = err.ctr,
        scale = err.scl
      )[, 1]
    ) # store scaled variable

  # densities by individual pattern
  epidens.largep <- mclapply(
    largep.list,
    \(x) {
      getepidens(
        x,
        data_fixed,
        response = response,
        runquiet = T,
        remove.tails = F,
        probs = c(0.1, 0.9),
        # continuousDiff = T,
        integrateToOne = F,
        # monotone = T,
        M = 5,
        unimodal = F
      )
    },
    mc.cores = mc
  )
  # store status
  # status.rep <- sapply(epidens.largep,
  # \(x) (paste("status:",x$status,"iterations:",x$iterations,", message: ",x$message,"\n")))

  # Markov clustering
  # Distance matrix
  # hypodist <- proxy::dist(epidens.largep, method = distance_epi,
  #                         diag = T,
  #                         upper = T) # runs too slow, explore alternatives

  hypodist <- dist.custom.mc(epidens.largep, distance_epi, mc.cores = mc)

  # prepare matrix for MCL
  adjacency <- clust_weights(hypodist, ...)
  # run MC with function
  pattern_mcl <- mcl(as.matrix(adjacency), addLoops = T)

  # add cluster for large patterns
  pattern.size$cluster <- NA
  pattern.size$cluster[pattern.size$size == "large"] <- pattern_mcl$Cluster

  # join small clusters
  pattern.small <- pattern.size %>%
    select(-cluster) %>%
    filter(size == "small") %>%
    stringdist_left_join(
      # fuzzy join with distance in characters
      x = .,
      y = pattern.size %>% filter(size == "large") %>% select(pattern, cluster),
      by = c("pattern"),
      max_dist = 1
    )

  # if several matches are found, take the most popular
  pattern.small <- pattern.small %>%
    group_by(pattern.x, cluster) %>%
    summarise(n = n()) %>%
    arrange(pattern.x, desc(n)) %>%
    summarise(cluster = first(cluster), .groups = "drop")

  # final clustering of patterns
  pattern.size.complete <- pattern.size %>%
    left_join(pattern.small, by = c("pattern" = "pattern.x")) %>%
    mutate(
      cluster = ifelse(is.na(cluster.x), cluster.y, cluster.x) %>%
        factor() %>%
        as.numeric()
    )

  pattern.size.complete
}


# get_skeleton_pts <- function(
#     wind.bpa.scen,
#     lead,
#     probs){
#   # filter data
#   wind.filtered <- wind.bpa.scen %>%
#     mutate(hour.lead=1:nrow(wind.bpa.scen)) %>%
#     filter(hour.lead == lead)
#   regime.ind <- regimelabel(wind.filtered$forecast,alpha = 0.4, data = wind.bpa2$forecast)
#   # debug(regimelabel)
#   wind.subset <- wind.bpa2 %>%
#     slice(regime.ind) %>%
#     filter(cluster == wind.filtered$cluster)
#
#   # estimate density
#   scaled.err <- wind.subset$err.cf %>% scale()
#   err.ctr <- attr(scaled.err,which = "scaled:center")
#   err.scl <- attr(scaled.err,which = "scaled:scale")
#
#   epidens.est <- expepi(scaled.err, softinfo = setup.softinfo(N = 10))
#
#   x <- epidens.est$x.pts
#   y <- epidens.est$y.est
#   nn <- length(x)
#   Fx <- cumsum(y * c(0, diff(x)))
#   Fx <- Fx/Fx[nn]
#
#   x <- epidens.est$x.pts * err.scl + err.ctr
#   y <- epidens.est$y.est / err.scl
#
#   dens.est <- density(wind.subset$err.cf)
#   dens.est$x <- x
#   dens.est$y <- y
#   quantile(dens.est, probs = probs)
#
#
#   # truncate if necessary
#   lwr.bound <- -wind.filtered$forecast
#   if (sum(x <= lwr.bound) > 0){
#     xprime <- c(lwr.bound-1,lwr.bound,x[x > lwr.bound])
#     discrete.prob <- Fx[max(which(x< lwr.bound))]
#     # yprime <- c(discrete.prob, y[x > lwr.bound]*(1-discrete.prob))
#     yprime <- c(0,discrete.prob, y[x > lwr.bound])
#   } else{
#     xprime <- x
#     discrete.prob <- 0
#     yprime <- y
#   }
#
#
#   # skeleton points
#   dens.est$x <- xprime
#   dens.est$y <- yprime
#
#   quants0 <- quantile(dens.est, probs = probs) %>%
#     pmax(lwr.bound)
#
#   skl.len <- length(quants0)-1
#   skeleton <-
#     sapply(1:skl.len,
#            \(i) {
#              skel.ind <- between(xprime,quants0[i],quants0[i+1])
#              part.mean <- sum( xprime[skel.ind] * yprime[skel.ind])/sum(yprime[skel.ind])
#              part.mean
#            })
#   skeleton
#
#   result = list(skeleton.pts = skeleton, epidens.est = epidens.est)
# }

get_skeleton_pts <- function(
  t1,
  full.data,
  response = "forecast.cf_orig",
  h,
  # day_part_separators,
  probs,
  probs.rel = FALSE
) {
  # look ahead window
  # h <- max(day_part_separators)

  # filter data
  wind.filtered <- full.data %>%
    filter(time >= t1, time <= t1 + hours(h)) %>%
    mutate(hour.lead = 0:h) %>%
    filter(hour.lead == h)

  regime.ind <- regimelabel(
    wind.filtered[[response]],
    alpha = 0.4,
    data = full.data[[response]]
  )
  # debug(regimelabel)

  wind.subset <- full.data %>%
    slice(regime.ind) %>%
    filter(cluster == wind.filtered$cluster)

  # estimate density
  scaled.err <- wind.subset$err.cf %>% scale()
  err.ctr <- attr(scaled.err, which = "scaled:center")
  err.scl <- attr(scaled.err, which = "scaled:scale")

  epidens.est <- expepi(
    scaled.err,
    softinfo = setup.softinfo(N = 10),
    postproc.controls = postproc.control(pic.types = NULL)
  )

  x <- epidens.est$x.pts
  y <- epidens.est$y.est
  nn <- length(x)
  Fx <- cumsum(y * c(0, diff(x)))
  Fx <- Fx / Fx[nn]

  x <- epidens.est$x.pts * err.scl + err.ctr
  y <- epidens.est$y.est / err.scl

  dens.est <- density(wind.subset$err.cf)
  dens.est$x <- x
  dens.est$y <- y
  # quantile(dens.est, probs = probs)

  # browser()
  # truncate if necessary
  lwr.bound <- -wind.filtered[[response]]
  upr.bound <- 1 - wind.filtered[[response]]

  if (sum(x <= lwr.bound) > 0) {
    # xprime <- c(-1,lwr.bound,x[x > lwr.bound])
    # discrete.prob <- Fx[max(which(x< lwr.bound))]
    # yprime <- c(0,discrete.prob, y[x > lwr.bound])

    xprime <- c(lwr.bound, x[x > lwr.bound])
    discrete.prob <- Fx[max(which(x < lwr.bound))]
    yprime <- c(discrete.prob, y[x > lwr.bound])
  } else {
    xprime <- x
    discrete.prob <- 0
    yprime <- y
  }

  if (sum(x >= upr.bound) > 0) {
    xprime2 <- c(xprime[xprime < upr.bound], upr.bound, 1)
    discrete.prob <- 1 - Fx[max(which(xprime > upr.bound))]
    yprime2 <- c(yprime[xprime < upr.bound], discrete.prob, 0)
  } else {
    xprime2 <- xprime
    discrete.prob <- 0
    yprime2 <- yprime
  }

  # skeleton points
  dens.est$x <- xprime2
  dens.est$y <- yprime2

  # browser()
  quants0 <- quantile(dens.est, probs = probs) %>%
    pmax(lwr.bound)

  skl.len <- length(quants0) - 1
  skeleton <-
    sapply(1:skl.len, \(i) {
      skel.ind <- between(xprime2, quants0[i], quants0[i + 1])
      part.mean <- sum(xprime2[skel.ind] * yprime2[skel.ind]) /
        sum(yprime2[skel.ind])
      part.mean
    })
  # skeleton

  result <- list(
    skeleton.pts = skeleton,
    quantiles = quants0,
    # quant.rel = quant.rel,
    epidens.est = epidens.est
  )

  if (probs.rel) {
    # 2nd set of quantiles
    probs.rel = seq(0.05, 0.95, by = 0.05) # probs for reliability diagram
    alpha <- 1 - probs.rel # alpha level
    quant_vector <- c(alpha / 2, 1 - alpha / 2) # quantiles for symmetric confidence interval
    quant.rel <- quantile(dens.est, probs = quant_vector)
    result$quant.rel <- quant.rel
  }

  result
}

# readwindyr <- function(year,
#                        path = "../data/",
#                        skip.lines = 23,
#                        column_names = c("time","forecast","actuals",
#                                         "load","hydro","thermal","netinter")){
#   file = paste(path,"WindGenTotalLoadYTD_",year,".xls",sep = "")
#   table <- read_xls(file, skip = skip.lines,
#                     col_names = column_names) %>% # read first semester
#     bind_rows(
#       read_xls(file, skip = skip.lines, sheet = 2,
#                col_names = column_names) # read second semester
#     ) %>%
#     mutate(time=as.POSIXct(time, tz = "PST", format = "%m/%d/%y %H:%M")) %>% # time format
#     group_by(time) %>%
#     summarise(forecast = first(forecast),actuals = first(actuals)) # remove extra data due to daylight savings time
#   table
# }

base_q_plot <- function(
  plot.dat,
  print.plot = TRUE
) {
  # browser()
  p <- plot.dat %>%
    ggplot() +
    geom_line(aes(time, value, col = name)) +
    theme_bw() +
    labs(col = "") +
    geom_point(
      aes(time, scenario.point),
      data = plot.dat #%>% na.omit()
    ) +
    ylab("Wind Power Output %") +
    coord_cartesian(ylim = c(0, 1)) +
    theme(
      legend.position = "bottom",
      # legend.position = "inside",legend.position.inside = c(.3,.2),
      legend.box = "horizontal",
      legend.background = element_rect(fill = NA, color = NA),
      legend.key = element_rect(fill = NA, color = NA)
    ) +
    scale_color_aaas()

  if (print.plot) {
    print(p)
  }
  invisible(p)
}

nonpar_scen_1day <- function(
  t1,
  day_part_separators = seq(0, 23),
  probs,
  data, # must have patterns / clusters
  plot.quant = TRUE,
  save.qplot = FALSE,
  fig.path = "~/Documents/proj2/fig_nonpar",
  sample.path = "~/Documents/proj2/sample",
  mod.code = "r_err_f_nonpar_t_dr_t",
  probs.rel = TRUE,
  save.samples = TRUE,
  parallel.run = TRUE, # run in parallel
  mc.cores = detectCores() - 2, # Default value for mc.cores
  ...
) {
  #

  # t1 <- as.POSIXct("2013-06-20 0:01")
  # probs <- seq(0.0, 1, by = 0.05)
  # day_part_separators <- seq(1,24,2)
  # day_part_separators <- c(1,12,24)

  # get forecast at DPS
  wind.bpa.scen <- data %>%
    filter(
      time >= t1 + hours(min(day_part_separators)),
      time <= t1 + hours(max(day_part_separators))
    ) %>%
    mutate(hour.lead = min(day_part_separators):max(day_part_separators))

  # find skelenton points
  skel.info <- if (parallel.run) {
    parallel::mclapply(
      day_part_separators,
      \(lead) {
        get_skeleton_pts(
          t1 = t1,
          full.data = data,
          h = lead,
          probs = probs,
          probs.rel = probs.rel
        )
      },
      mc.cores = mc.cores
    )
  } else {
    lapply(
      day_part_separators,
      \(lead) {
        get_skeleton_pts(
          t1 = t1,
          full.data = data,
          h = lead,
          probs = probs,
          probs.rel = probs.rel
        )
      }
    )
  }

  # Generate short names for the skeleton points
  short_names <- paste0(
    "p_",
    sprintf("%02.0f", prob_cut_points[-length(prob_cut_points)] * 100),
    "_",
    sprintf("%02.0f", prob_cut_points[-1] * 100)
  )

  # format skeleton points
  skel.dat <- lapply(skel.info, \(x) x$skeleton.pts) %>%
    unlist() %>%
    matrix(ncol = length(day_part_separators)) %>%
    t() %>%
    as.data.frame() %>%
    setNames(., short_names) %>%
    mutate(hour.lead = day_part_separators) %>%
    pivot_longer(
      cols = matches("p_"),
      names_to = "scen",
      values_to = "skel"
    ) %>%
    mutate(
      name = "forecast.cf_orig",
      scen = factor(
        scen,
        levels = short_names,
        labels = short_names
      )
    )

  # browser()
  # format quantiles
  quant.dat <- lapply(
    skel.info,
    \(x) x$quant.rel
  ) %>%
    bind_rows() %>%
    # unlist() %>% matrix(ncol =length(day_part_separators)) %>% t() %>%
    # as.data.frame() %>%
    setNames(., paste0("quantile_", names(.))) %>%
    mutate(hour.lead = day_part_separators) %>%
    # pivot_longer(cols = matches("quantile_"), names_to = "scen",
    #              values_to = "quant") %>%
    mutate(
      name = "forecast.cf_orig",
      scen = factor(
        scen
      )
    )

  # skeleton plot data
  skel.plot.data <- wind.bpa.scen %>%
    pivot_longer(cols = c(forecast.cf_orig, actuals.cf)) %>%
    # filter(hour.lead %in% day_part_separators) %>%
    left_join(skel.dat, by = c("hour.lead", "name")) %>%
    mutate(scenario.point = skel + value)

  skel.data.samples <- skel.plot.data %>%
    select(time, scen, scenario.point) %>%
    mutate(time = paste0("t_", format(time, "%Y-%m-%d_%H"))) %>%
    group_by(scen) %>%
    pivot_wider(names_from = time, values_from = scenario.point) %>%
    na.omit()

  # quantiles plot data
  quant.plot.data <- wind.bpa.scen %>%
    # pivot_longer(cols = c(forecast.cf_orig,actuals.cf)) %>%
    # filter(hour.lead %in% day_part_separators) %>%
    # left_join(quant.dat, by = c("hour.lead","name")) %>%
    left_join(quant.dat, by = c("hour.lead")) %>%
    # mutate(scenario.point = quant + value)
    mutate(across(matches("quantile_"), ~ . + forecast.cf_orig))

  # browser()
  p.q <- base_q_plot(skel.plot.data, print.plot = plot.quant)

  if (save.qplot) {
    ggsave(
      file.path(fig.path, paste0(mod.code, "_skel.png")),
      p.q,
      width = 5,
      height = 3.5
    )
  }

  if (save.samples) {
    # save skeleton points scenarios in same format as inla samples
    skel.data.samples %>%
      write.csv(
        file = file.path(
          sample.path,
          paste0(mod.code, t1, ".csv")
        )
      )

    quant.plot.data %>%
      write.csv(
        file = file.path(
          sample.path,
          paste0(mod.code, "quant_t", t1, ".csv")
        )
      )
  } else {
    result <- list(
      skel.dat = skel.plot.data,
      quant.dat = quant.plot.data
    )
    invisible(result)
  }
}


nonpar_scen_ndays <- function(
  inidate,
  n.days,
  day_part_separators,
  prob_cut_points
) {
  # browser()
  # day sequence to predict day-ahead power
  time_seq <- seq(
    from = inidate,
    to = inidate + n.days * 24 * 3600,
    by = "day"
  )

  # record initial time
  start_time <- Sys.time()

  # debug(nonpar_scen_1day)
  lapply(
    time_seq,
    \(time) {
      nonpar_scen_1day(
        t1 = time,
        day_part_separators = day_part_separators,
        probs = prob_cut_points,
        data = data.scaled, # must have patterns / clusters
        plot.quant = FALSE,
        save.qplot = FALSE,
        parallel.run = TRUE
      )
      # print(paste("Day", time, "completed"))
    }
  )

  end_time <- Sys.time()
  run_time <- end_time - start_time
  print(run_time)
}

nonpar.scoring.reliability <- function(
  model.name,
  time_seq,
  sample.path = file.path("~/Documents/proj2", "sample"),
  ...
) {
  # get response from model code
  response = "actuals.cf"

  probs = seq(0.05, 0.95, by = 0.05) # probs for reliability diagram
  alpha <- 1 - probs # alpha level
  quant_vector <- c(alpha / 2, 1 - alpha / 2) # quantiles for symmetric confidence interval

  # get quantiles
  q.reliability <- lapply(
    time_seq,
    \(day) {
      fread(
        file.path(
          sample.path,
          paste0(model.name, "quant_t", format(day, "%Y-%m-%d"), ".csv")
        )
      )
    }
  ) %>%
    bind_rows()

  # get skeleton related scores
  # skel.scores <- lapply(
  #   time_seq,
  #   \(day) {
  #     nonpar_scen_scoring(
  #       sample.path = sample.path,
  #       model.name = model.name,
  #       day,
  #       response
  #     )
  #   }) %>%
  #   bind_rows()

  skel.scores <- lapply(
    time_seq,
    \(day) {
      tryCatch(
        {
          nonpar_scen_scoring(
            sample.path = sample.path,
            model.name = model.name,
            day,
            response
          )
        },
        error = function(e) {
          message("Error encountered on day: ", day)
          print(e) # Print error message
          browser() # Enter debug mode
          # NULL      # Return NULL or other placeholder if needed
          nonpar_scen_scoring(
            sample.path = sample.path,
            model.name = model.name,
            day,
            response
          )
        }
      )
    }
  ) %>%
    bind_rows()

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

  # calculate all coverages by day
  daily.scores.model <- q.reliability %>%
    left_join(
      skel.scores %>% select(time, any_of(stats.vec)),
      by = "time"
    ) %>%
    group_by(date) %>%
    summarise(
      rmse = sqrt(mean((.data[[response]] - mean)^2)),
      across(!!stats.vec, mean),
      coverage05 = mean(between(
        .data[[response]],
        .data[[quant_low[1]]],
        .data[[quant_high[1]]]
      )),
      coverage10 = mean(between(
        .data[[response]],
        .data[[quant_low[2]]],
        .data[[quant_high[2]]]
      )),
      coverage15 = mean(between(
        .data[[response]],
        .data[[quant_low[3]]],
        .data[[quant_high[3]]]
      )),
      coverage20 = mean(between(
        .data[[response]],
        .data[[quant_low[4]]],
        .data[[quant_high[4]]]
      )),
      coverage25 = mean(between(
        .data[[response]],
        .data[[quant_low[5]]],
        .data[[quant_high[5]]]
      )),
      coverage30 = mean(between(
        .data[[response]],
        .data[[quant_low[6]]],
        .data[[quant_high[6]]]
      )),
      coverage35 = mean(between(
        .data[[response]],
        .data[[quant_low[7]]],
        .data[[quant_high[7]]]
      )),
      coverage40 = mean(between(
        .data[[response]],
        .data[[quant_low[8]]],
        .data[[quant_high[8]]]
      )),
      coverage45 = mean(between(
        .data[[response]],
        .data[[quant_low[9]]],
        .data[[quant_high[9]]]
      )),
      coverage50 = mean(between(
        .data[[response]],
        .data[[quant_low[10]]],
        .data[[quant_high[10]]]
      )),
      coverage55 = mean(between(
        .data[[response]],
        .data[[quant_low[11]]],
        .data[[quant_high[11]]]
      )),
      coverage60 = mean(between(
        .data[[response]],
        .data[[quant_low[12]]],
        .data[[quant_high[12]]]
      )),
      coverage65 = mean(between(
        .data[[response]],
        .data[[quant_low[13]]],
        .data[[quant_high[13]]]
      )),
      coverage70 = mean(between(
        .data[[response]],
        .data[[quant_low[14]]],
        .data[[quant_high[14]]]
      )),
      coverage75 = mean(between(
        .data[[response]],
        .data[[quant_low[15]]],
        .data[[quant_high[15]]]
      )),
      coverage80 = mean(between(
        .data[[response]],
        .data[[quant_low[16]]],
        .data[[quant_high[16]]]
      )),
      coverage85 = mean(between(
        .data[[response]],
        .data[[quant_low[17]]],
        .data[[quant_high[17]]]
      )),
      coverage90 = mean(between(
        .data[[response]],
        .data[[quant_low[18]]],
        .data[[quant_high[18]]]
      )),
      coverage95 = mean(between(
        .data[[response]],
        .data[[quant_low[19]]],
        .data[[quant_high[19]]]
      )) #,
      # coverage005 = mean(.data[[response]]<`quantile_0.5%`)
    ) %>%
    mutate(
      model = model.name,
      response = response
    )

  # same calculation by hour
  hour.scores.model <- q.reliability %>%
    left_join(
      skel.scores %>% select(time, any_of(stats.vec)),
      by = "time"
    ) %>%
    group_by(hour) %>%
    summarise(
      rmse = sqrt(mean((.data[[response]] - mean)^2)),
      across(!!stats.vec, mean),
      coverage05 = mean(between(
        .data[[response]],
        .data[[quant_low[1]]],
        .data[[quant_high[1]]]
      )),
      coverage10 = mean(between(
        .data[[response]],
        .data[[quant_low[2]]],
        .data[[quant_high[2]]]
      )),
      coverage15 = mean(between(
        .data[[response]],
        .data[[quant_low[3]]],
        .data[[quant_high[3]]]
      )),
      coverage20 = mean(between(
        .data[[response]],
        .data[[quant_low[4]]],
        .data[[quant_high[4]]]
      )),
      coverage25 = mean(between(
        .data[[response]],
        .data[[quant_low[5]]],
        .data[[quant_high[5]]]
      )),
      coverage30 = mean(between(
        .data[[response]],
        .data[[quant_low[6]]],
        .data[[quant_high[6]]]
      )),
      coverage35 = mean(between(
        .data[[response]],
        .data[[quant_low[7]]],
        .data[[quant_high[7]]]
      )),
      coverage40 = mean(between(
        .data[[response]],
        .data[[quant_low[8]]],
        .data[[quant_high[8]]]
      )),
      coverage45 = mean(between(
        .data[[response]],
        .data[[quant_low[9]]],
        .data[[quant_high[9]]]
      )),
      coverage50 = mean(between(
        .data[[response]],
        .data[[quant_low[10]]],
        .data[[quant_high[10]]]
      )),
      coverage55 = mean(between(
        .data[[response]],
        .data[[quant_low[11]]],
        .data[[quant_high[11]]]
      )),
      coverage60 = mean(between(
        .data[[response]],
        .data[[quant_low[12]]],
        .data[[quant_high[12]]]
      )),
      coverage65 = mean(between(
        .data[[response]],
        .data[[quant_low[13]]],
        .data[[quant_high[13]]]
      )),
      coverage70 = mean(between(
        .data[[response]],
        .data[[quant_low[14]]],
        .data[[quant_high[14]]]
      )),
      coverage75 = mean(between(
        .data[[response]],
        .data[[quant_low[15]]],
        .data[[quant_high[15]]]
      )),
      coverage80 = mean(between(
        .data[[response]],
        .data[[quant_low[16]]],
        .data[[quant_high[16]]]
      )),
      coverage85 = mean(between(
        .data[[response]],
        .data[[quant_low[17]]],
        .data[[quant_high[17]]]
      )),
      coverage90 = mean(between(
        .data[[response]],
        .data[[quant_low[18]]],
        .data[[quant_high[18]]]
      )),
      coverage95 = mean(between(
        .data[[response]],
        .data[[quant_low[19]]],
        .data[[quant_high[19]]]
      )) #,
      # coverage005 = mean(.data[[response]]<`quantile_0.5%`)
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
    hour = hour.scores.model
  ))
}


nonpar_scen_scoring <- function(
  sample.path = file.path("~/Documents/proj2", "sample"),
  model.name,
  t1,
  response,
  probs = c(0.025, 0.5, 0.975)
) {
  scen.tbl <- file.path(
    sample.path,
    paste0(model.name, date(t1), ".csv")
  ) %>%
    fread(.) %>%
    as.data.frame() %>%
    select(starts_with("t_"))

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
  function_list <- list(
    mean = mean,
    sd = sd,
    quantile = quantile
  )

  # Call the function with desired quantiles
  # result <- functions_to_scenarios(sim_data, function_list, probs = c(0.05, 0.5, 0.95))
  result <- functions_to_scenarios(
    scen.tbl,
    function_list,
    probs = probs
  ) %>%
    mutate(
      # Remove "t_" prefix and convert to POSIXct, handling NA text entries
      time = if_else(
        str_detect(time, "^NA$"), # Check if the value is "NA"
        NA_POSIXct_, # Keep it as NA if it's "NA"
        as.POSIXct(str_remove(time, "^t_"), format = "%Y-%m-%d_%H", tz = "PST")
      )
    )
  # browser()
  crps <- scoringRules::crps_sample(
    y = subsample %>% pull(!!response), # day-ahead response
    dat = scen.tbl %>%
      select(1:nrow(subsample)) %>% # select 24-hours
      as.matrix() %>%
      t() # transpose so that crps works
  )

  energy.s <- scoringRules::es_sample(
    y = subsample %>% pull(!!response), # day-ahead response
    dat = scen.tbl %>%
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
              pull(!!response), # day-ahead response
            dat = scen.tbl %>%
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
