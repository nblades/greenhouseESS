############################
# Pooled variance/covariance traces & plotting
############################

# pooled variance traces

s2p_trace_one <- function(x, Tmax = NULL, progress = TRUE) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  N <- length(x)
  if (is.null(Tmax)) Tmax <- floor(N / 10)
  if (Tmax < 1) stop("Not enough data: N must be >= 10.")
  
  Tvec <- seq_len(Tmax)
  np   <- N / Tvec                      # continuous, not floor()
  s2p  <- numeric(Tmax)
  
  if (progress) {
    message(sprintf("[variance] N=%d, Tmax=%d", N, Tmax))
    flush.console()
    pb <- txtProgressBar(min = 1, max = Tmax, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  
  # T = 1
  s2p[1] <- var(x)
  if (progress) setTxtProgressBar(pb, 1L)
  
  if (Tmax >= 2) {
    for (T in 2:Tmax) {
      s2_t <- numeric(T)
      for (i in seq_len(T)) {
        idx <- seq.int(i, N, by = T)
        s2_t[i] <- if (length(idx) >= 2) var(x[idx]) else NA_real_
      }
      s2p[T] <- mean(s2_t, na.rm = TRUE)
      if (progress) setTxtProgressBar(pb, T)
    }
  }
  
  tibble(T = Tvec, np = np, s2p = s2p)
}

compute_s2p_traces <- function(data,
                               vars = c("co2_ppm", "humidity_pct", "soil_temp_F"),
                               progress = TRUE) {
  stopifnot(all(vars %in% names(data)))
  if (progress) {
    message(sprintf("Computing pooled variances for %d variables ...", length(vars)))
    flush.console()
  }
  
  out <- lapply(seq_along(vars), function(j) {
    v <- vars[j]
    if (progress) {
      message(sprintf("  [%d/%d] %s", j, length(vars), v))
      flush.console()
    }
    ans <- s2p_trace_one(data[[v]], progress = progress)
    ans$variable <- v
    ans
  })
  
  bind_rows(out) |>
    mutate(variable = factor(.data$variable, levels = vars))
}

# pooled covariance traces (with progress) 

compute_pairwise_cov <- function(data,
                                 vars = c("co2_ppm", "humidity_pct", "soil_temp_F"),
                                 progress = TRUE) {
  stopifnot(all(vars %in% names(data)))
  pairs <- combn(vars, 2, simplify = FALSE)
  
  out <- lapply(seq_along(pairs), function(k) {
    pr <- pairs[[k]]
    a  <- pr[1]
    b  <- pr[2]
    if (progress) message(sprintf("  [%d/%d] %s × %s", k, length(pairs), a, b))
    
    x <- as.numeric(data[[a]])
    y <- as.numeric(data[[b]])
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    n <- length(x)
    if (n < 2) return(NULL)
    
    Tmax <- floor(n / 10)
    Tvec <- seq_len(Tmax)
    np   <- n / Tvec
    sXYp <- numeric(Tmax)
    
    if (progress) {
      pb <- txtProgressBar(min = 1, max = Tmax, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    
    # T = 1
    sXYp[1] <- cov(x, y)
    if (progress) setTxtProgressBar(pb, 1L)
    
    if (Tmax >= 2) {
      for (T in 2:Tmax) {
        cov_t <- numeric(T)
        for (i in seq_len(T)) {
          idx <- seq.int(i, n, by = T)
          cov_t[i] <- if (length(idx) >= 2) cov(x[idx], y[idx]) else NA_real_
        }
        sXYp[T] <- mean(cov_t, na.rm = TRUE)
        if (progress) setTxtProgressBar(pb, T)
      }
    }
    
    tibble(
      var_x    = a,
      var_y    = b,
      pair_key = paste(a, b, sep = "_"),
      T        = Tvec,
      np       = np,
      sXYp     = sXYp
    )
  })
  
  bind_rows(out)
}

# labels & T-masking

pretty_var <- function(v) {
  c(
    co2_ppm      = "CO[2]",
    humidity_pct = "Humidity",
    soil_temp_F  = "Soil~Temperature"
  )[v]
}

# stacked, two-line facet label for covariance: atop(top, bottom)
pair_label_stacked <- function(x, y) {
  paste0("atop(", pretty_var(x), ", ", pretty_var(y), ")")
}

# Mask a y-column based on T windows; everything else untouched.
mask_T_windows <- function(df,
                           T_col,
                           y_col,
                           highlight_Ts,
                           highlight_T_halfwidth) {
  # If nothing to mask, just copy through
  if (is.null(highlight_Ts) || length(highlight_Ts) == 0) {
    df[[paste0(y_col, "_mask")]] <- df[[y_col]]
    return(df)
  }
  
  # Recycle halfwidth if scalar
  if (length(highlight_T_halfwidth) == 1L &&
      length(highlight_Ts) > 1L) {
    highlight_T_halfwidth <- rep(highlight_T_halfwidth, length(highlight_Ts))
  }
  stopifnot(length(highlight_Ts) == length(highlight_T_halfwidth))
  
  Tvals <- df[[T_col]]
  bad   <- rep(FALSE, length(Tvals))
  
  for (k in seq_along(highlight_Ts)) {
    T0 <- highlight_Ts[k]
    hw <- highlight_T_halfwidth[k]
    lo <- T0 - hw
    hi <- T0 + hw
    bad <- bad | (Tvals >= lo & Tvals <= hi)
  }
  
  y <- df[[y_col]]
  df[[paste0(y_col, "_mask")]] <- ifelse(bad, NA_real_, y)
  df
}

# pooled variance (left) + covariance (right), T-masked 

plot_pooled_var_cov <- function(
    var_traces,
    cov_traces,
    col = if (exists("gh_cols")) gh_cols("moody") else "#36648B",
    filename = NULL,
    width = 10, height = 8, dpi = 300,
    x_limits = c(10, 300),
    x_breaks  = c(10, 50, 100, 150, 200, 250, 300),
    
    # mask windows in T (minutes)
    # 8, 12, 16, 24, 36, 48, 60, 72 hours
    highlight_Ts          = c(480, 720, 960, 1440, 2160, 2880, 3600, 4320),
    highlight_T_halfwidth = c(20,  25,  20,   95,   40,   175,  60,   60),
    
    # show selected T's as diamonds (from choose_greenhouse_Ts())
    res_T          = NULL,          # list(var_choices, cov_choices)
    sel_method     = "roll",        # which method's T to show
    sel_point_col  = if (exists("gh_cols")) gh_cols("orange") else "#228B22"
) {
  stopifnot(all(c("T", "np", "s2p",  "variable") %in% names(var_traces)))
  stopifnot(all(c("T", "np", "sXYp", "var_x", "var_y", "pair_key") %in% names(cov_traces)))
  
  # Filter to the n-range you care about, then mask in T-space
  vdat <- var_traces |>
    filter(np > x_limits[1], np < x_limits[2]) |>
    mask_T_windows(
      T_col = "T",
      y_col = "s2p",
      highlight_Ts          = highlight_Ts,
      highlight_T_halfwidth = highlight_T_halfwidth
    )
  
  cdat <- cov_traces |>
    mutate(pair_label = pair_label_stacked(var_x, var_y)) |>
    filter(np > x_limits[1], np < x_limits[2]) |>
    mask_T_windows(
      T_col = "T",
      y_col = "sXYp",
      highlight_Ts          = highlight_Ts,
      highlight_T_halfwidth = highlight_T_halfwidth
    )
  
  # T = 1 references (unmasked)
  t1_var <- var_traces |>
    filter(T == 1) |>
    select(variable, s2p)
  
  t1_cov <- cov_traces |>
    filter(T == 1) |>
    mutate(pair_label = pair_label_stacked(var_x, var_y)) |>
    select(pair_label, sXYp)
  
  # Selected Ts (ESS choices) if supplied
  var_pts <- NULL
  cov_pts <- NULL
  
  if (!is.null(res_T)) {
    if ("var_choices" %in% names(res_T)) {
      var_pts <- res_T$var_choices |>
        filter(method == sel_method,
               !is.na(np_selected),
               !is.na(estimate)) |>
        select(variable, np_selected, estimate) |>
        distinct() |>
        filter(np_selected > x_limits[1], np_selected < x_limits[2])
      
      if (!is.null(vdat$variable)) {
        var_pts$variable <- factor(var_pts$variable, levels = levels(vdat$variable))
      }
    }
    
    if ("cov_choices" %in% names(res_T)) {
      cov_pts <- res_T$cov_choices |>
        filter(method == sel_method,
               !is.na(np_selected),
               !is.na(estimate)) |>
        mutate(pair_label = pair_label_stacked(var_x, var_y)) |>
        select(pair_label, np_selected, estimate) |>
        distinct() |>
        filter(np_selected > x_limits[1], np_selected < x_limits[2])
    }
  }
  
  var_labeller <- as_labeller(
    c(
      co2_ppm      = "CO[2]~'(ppm)'",
      humidity_pct = "'Humidity (%)'",
      soil_temp_F  = "'Soil Temp ('*degree*'F)'"
    ),
    label_parsed
  )
  
  base_theme <- theme_bw() +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_blank(),
      strip.placement   = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "plain"),
      plot.title        = element_text(hjust = 0.5, face = "bold")
    )
  
  # left: pooled variance (masked)
  p_var <- ggplot(vdat, aes(x = np, y = s2p_mask)) +
    geom_line(color = col, na.rm = TRUE) +
    {
      if (!is.null(var_pts) && nrow(var_pts) > 0)
        geom_point(
          data = var_pts,
          inherit.aes = FALSE,
          aes(x = np_selected, y = estimate),
          shape  = 23,
          size   = 2.8,
          stroke = 0.4,
          color  = "black",
          fill   = sel_point_col
        )
    } +
    geom_hline(
      data = t1_var,
      aes(yintercept = s2p),
      inherit.aes = FALSE,
      color = "#63B8FF",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ variable,
      ncol  = 1,
      scales = "free_y",
      strip.position = "left",
      labeller       = var_labeller
    ) +
    scale_x_sqrt(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_breaks,
      name   = expression(n[t])
    ) +
    labs(y = "", title = expression("Pooled Variance " * s[i*","*"ESS"]^2)) +
    base_theme
  
  # right: pooled covariance (masked)
  p_cov <- ggplot(cdat, aes(x = np, y = sXYp_mask)) +
    geom_line(color = col, na.rm = TRUE) +
    {
      if (!is.null(cov_pts) && nrow(cov_pts) > 0)
        geom_point(
          data = cov_pts,
          inherit.aes = FALSE,
          aes(x = np_selected, y = estimate),
          shape  = 23,
          size   = 2.8,
          stroke = 0.4,
          color  = "black",
          fill   = sel_point_col
        )
    } +
    geom_hline(
      data = t1_cov,
      aes(yintercept = sXYp),
      inherit.aes = FALSE,
      color = "#63B8FF",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ pair_label,
      ncol  = 1,
      scales = "free_y",
      strip.position = "left",
      labeller       = label_parsed
    ) +
    scale_x_sqrt(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_breaks,
      name   = expression(n[t])
    ) +
    labs(y = "", title = expression("Pooled Covariance " * s[i*j*","*"ESS"])) +
    base_theme
  
  combo <- (p_var | p_cov) + plot_layout(widths = c(1, 1))
  
  if (!is.null(filename)) {
    ggsave(filename, combo, width = width, height = height, dpi = dpi)
  }
  
  combo
}

# single-trace helpers (no masking) 

plot_single_var_trace <- function(var_traces,
                                  variable,
                                  col = if (exists("gh_cols")) gh_cols("moody") else "#36648B",
                                  x_limits = c(10, 300),
                                  x_breaks = c(10, 50, 100, 150, 200, 250, 300)) {
  dat <- var_traces |>
    filter(variable == !!variable,
           np > x_limits[1],
           np < x_limits[2])
  
  ggplot(dat, aes(x = np, y = s2p)) +
    geom_line(color = col) +
    scale_x_sqrt(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_breaks,
      name   = expression(n[t])
    ) +
    labs(
      y     = "",
      title = paste("Pooled variance trace:", variable)
    ) +
    theme_bw() +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

plot_single_cov_trace <- function(cov_traces,
                                  var_x,
                                  var_y,
                                  col = if (exists("gh_cols")) gh_cols("moody") else "#36648B",
                                  x_limits = c(10, 300),
                                  x_breaks = c(10, 50, 100, 150, 200, 250, 300)) {
  key1 <- paste(var_x, var_y, sep = "_")
  key2 <- paste(var_y, var_x, sep = "_")
  
  dat <- cov_traces |>
    filter(pair_key %in% c(key1, key2),
           np > x_limits[1],
           np < x_limits[2]) |>
    mutate(pair_label = pair_label_stacked(var_x, var_y))
  
  ggplot(dat, aes(x = np, y = sXYp)) +
    geom_line(color = col) +
    scale_x_sqrt(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_breaks,
      name   = expression(n[t])
    ) +
    labs(
      y     = "",
      title = paste("Pooled covariance trace:", var_x, "×", var_y)
    ) +
    theme_bw() +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

############################
# Plot for tuning T-mask windows (T-space only)
############################

plot_pooled_var_cov_T_mask_tuning <- function(
    var_traces,
    cov_traces,
    col = if (exists("gh_cols")) gh_cols("moody") else "#36648B",
    filename = NULL,
    width = 10, height = 8, dpi = 300,
    
    # x-axis in T
    T_limits = c(10, 5000),
    T_breaks = c(10, 100, 500, 1000, 2000, 4000),
    log10_x  = TRUE,
    
    # resonance T's (centers) and windows around them, in minutes
    highlight_Ts          = c(720, 1440, 2160, 2880, 3600, 4320),
    highlight_T_halfwidth = 100,
    
    shade_color = "hotpink",
    shade_alpha = 0.15
) {
  stopifnot(all(c("T", "np", "s2p",  "variable") %in% names(var_traces)))
  stopifnot(all(c("T", "np", "sXYp", "var_x", "var_y", "pair_key") %in% names(cov_traces)))
  
  vdat <- var_traces |>
    filter(T >= T_limits[1], T <= T_limits[2])
  
  cdat <- cov_traces |>
    filter(T >= T_limits[1], T <= T_limits[2]) |>
    mutate(pair_label = pair_label_stacked(var_x, var_y))
  
  t1_var <- var_traces |>
    filter(T == 1) |>
    select(variable, s2p)
  
  t1_cov <- cov_traces |>
    filter(T == 1) |>
    mutate(pair_label = pair_label_stacked(var_x, var_y)) |>
    select(pair_label, sXYp)
  
  if (length(highlight_Ts) > 0) {
    if (length(highlight_T_halfwidth) == 1L &&
        length(highlight_Ts) > 1L) {
      highlight_T_halfwidth <- rep(highlight_T_halfwidth, length(highlight_Ts))
    }
    stopifnot(length(highlight_Ts) == length(highlight_T_halfwidth))
    
    rects <- tibble(
      Tmin  = highlight_Ts - highlight_T_halfwidth,
      Tmax  = highlight_Ts + highlight_T_halfwidth,
      label = paste0("T≈", highlight_Ts)
    ) |>
      mutate(
        Tmin = pmax(Tmin, T_limits[1]),
        Tmax = pmin(Tmax, T_limits[2])
      ) |>
      filter(Tmin < Tmax)
  } else {
    rects <- tibble(Tmin = numeric(0), Tmax = numeric(0), label = character(0))
  }
  
  var_labeller <- as_labeller(
    c(
      co2_ppm      = "CO[2]~'(ppm)'",
      humidity_pct = "'Humidity (%)'",
      soil_temp_F  = "'Soil Temp ('*degree*'F)'"
    ),
    label_parsed
  )
  
  base_theme <- theme_bw() +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_blank(),
      strip.placement   = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "plain"),
      plot.title        = element_text(hjust = 0.5, face = "bold")
    )
  
  x_scale_fun <- if (log10_x) {
    scale_x_log10(
      limits = T_limits,
      breaks = T_breaks,
      labels = T_breaks,
      name   = "T"
    )
  } else {
    scale_x_continuous(
      limits = T_limits,
      breaks = T_breaks,
      labels = T_breaks,
      name   = "T"
    )
  }
  
  p_var <- ggplot(vdat, aes(x = T, y = s2p)) +
    {
      if (nrow(rects) > 0)
        geom_rect(
          data = rects,
          inherit.aes = FALSE,
          aes(xmin = Tmin, xmax = Tmax, ymin = -Inf, ymax = Inf),
          fill  = shade_color,
          alpha = shade_alpha
        )
    } +
    geom_line(color = col, na.rm = TRUE) +
    geom_hline(
      data = t1_var,
      aes(yintercept = s2p),
      inherit.aes = FALSE,
      color = "#63B8FF",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ variable,
      ncol = 1,
      scales = "free_y",
      strip.position = "left",
      labeller = var_labeller
    ) +
    x_scale_fun +
    labs(y = "", title = expression("Pooled Variance " * s[i*","*"ESS"]^2)) +
    base_theme
  
  p_cov <- ggplot(cdat, aes(x = T, y = sXYp)) +
    {
      if (nrow(rects) > 0)
        geom_rect(
          data = rects,
          inherit.aes = FALSE,
          aes(xmin = Tmin, xmax = Tmax, ymin = -Inf, ymax = Inf),
          fill  = shade_color,
          alpha = shade_alpha
        )
    } +
    geom_line(color = col, na.rm = TRUE) +
    geom_hline(
      data = t1_cov,
      aes(yintercept = sXYp),
      inherit.aes = FALSE,
      color = "#63B8FF",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ pair_label,
      ncol = 1,
      scales = "free_y",
      strip.position = "left",
      labeller = label_parsed
    ) +
    x_scale_fun +
    labs(y = "", title = expression("Pooled Covariance " * s[i*j*","*"ESS"])) +
    base_theme
  
  combo <- (p_var | p_cov) + plot_layout(widths = c(1, 1))
  
  if (!is.null(filename)) {
    ggsave(filename, combo, width = width, height = height, dpi = dpi)
  }
  
  combo
}

############################
# Pooled var/cov vs n_t with T-based shaded windows
############################

plot_pooled_var_cov_mask_tuning <- function(
    var_traces,
    cov_traces,
    col = if (exists("gh_cols")) gh_cols("moody") else "#36648B",
    filename = NULL,
    width = 10, height = 8, dpi = 300,
    x_limits = c(10, 300),
    x_breaks  = c(10, 50, 100, 150, 200, 250, 300),
    
    # T windows (minutes) that drive everything
    # e.g. 8, 12, 16, 24, 36, 48, 60, 72 hours
    highlight_Ts          = c(480, 720, 960, 1440, 2160, 2880, 3600, 4320),
    highlight_T_halfwidth = c(20,  25,  20,   95,   40,   175,  60,   60),
    
    shade_color = "hotpink",
    shade_alpha = 0.12,
    
    dense_ticks    = FALSE,
    n_dense_breaks = 12
) {
  stopifnot(all(c("T", "np", "s2p",  "variable") %in% names(var_traces)))
  stopifnot(all(c("T", "np", "sXYp", "var_x", "var_y", "pair_key") %in% names(cov_traces)))
  
  # keep full curves, just limit in n
  vdat <- var_traces |>
    dplyr::filter(np > x_limits[1], np < x_limits[2])
  
  cdat <- cov_traces |>
    dplyr::filter(np > x_limits[1], np < x_limits[2]) |>
    dplyr::mutate(pair_label = pair_label_stacked(var_x, var_y))
  
  # Estimate N from traces: np * T ≈ N
  N_est <- round(stats::median(var_traces$np * var_traces$T, na.rm = TRUE))
  
  # Build n_t rectangles from T windows only
  if (length(highlight_Ts) > 0) {
    if (length(highlight_T_halfwidth) == 1L &&
        length(highlight_Ts) > 1L) {
      highlight_T_halfwidth <- rep(highlight_T_halfwidth, length(highlight_Ts))
    }
    stopifnot(length(highlight_Ts) == length(highlight_T_halfwidth))
    
    rects <- purrr::map2_dfr(highlight_Ts, highlight_T_halfwidth, function(T0, hw) {
      T_low  <- max(1, T0 - hw)
      T_high <- max(T_low + 1e-8, T0 + hw)
      # Convert T-window to n_t window via N_est / T
      n_low  <- N_est / T_high
      n_high <- N_est / T_low
      
      tibble::tibble(
        xmin  = min(n_low, n_high),
        xmax  = max(n_low, n_high),
        label = paste0("T≈", T0)
      )
    }) |>
      dplyr::mutate(
        xmin = pmax(xmin, x_limits[1]),
        xmax = pmin(xmax, x_limits[2])
      ) |>
      dplyr::filter(xmin < xmax)
  } else {
    rects <- tibble::tibble(
      xmin  = numeric(0),
      xmax  = numeric(0),
      label = character(0)
    )
  }
  
  # T = 1 reference lines
  t1_var <- var_traces |>
    dplyr::filter(T == 1) |>
    dplyr::select(variable, s2p)
  
  t1_cov <- cov_traces |>
    dplyr::filter(T == 1) |>
    dplyr::mutate(pair_label = pair_label_stacked(var_x, var_y)) |>
    dplyr::select(pair_label, sXYp)
  
  # x-axis ticks
  x_breaks_use <- if (dense_ticks) {
    pretty(x_limits, n = n_dense_breaks)
  } else {
    x_breaks
  }
  
  # same labeller as main plot
  var_labeller <- as_labeller(
    c(
      co2_ppm      = "CO[2]~'(ppm)'",
      humidity_pct = "'Humidity (%)'",
      soil_temp_F  = "'Soil Temp ('*degree*'F)'"
    ),
    label_parsed
  )
  
  base_theme <- theme_bw() +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_blank(),
      strip.placement   = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "plain"),
      plot.title        = element_text(hjust = 0.5, face = "bold")
    )
  
  # left: pooled variance vs n_t, shaded windows (no masking)
  p_var <- ggplot(vdat, aes(x = np, y = s2p)) +
    {
      if (nrow(rects) > 0)
        geom_rect(
          data = rects,
          inherit.aes = FALSE,
          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
          fill  = shade_color,
          alpha = shade_alpha
        )
    } +
    geom_line(color = col, na.rm = TRUE) +
    geom_hline(
      data = t1_var,
      aes(yintercept = s2p),
      inherit.aes = FALSE,
      color = "#63B8FF",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ variable,
      ncol = 1,
      scales = "free_y",
      strip.position = "left",
      labeller = var_labeller
    ) +
    scale_x_sqrt(
      limits = x_limits,
      breaks = x_breaks_use,
      labels = x_breaks_use,
      name   = "n"        # or expression(n[t]) if you prefer
    ) +
    labs(y = "", title = expression("Pooled Variance " * s[i*","*"ESS"]^2)) +
    base_theme
  
  # right: pooled covariance vs n_t, shaded windows (no masking)
  p_cov <- ggplot(cdat, aes(x = np, y = sXYp)) +
    {
      if (nrow(rects) > 0)
        geom_rect(
          data = rects,
          inherit.aes = FALSE,
          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
          fill  = shade_color,
          alpha = shade_alpha
        )
    } +
    geom_line(color = col, na.rm = TRUE) +
    geom_hline(
      data = t1_cov,
      aes(yintercept = sXYp),
      inherit.aes = FALSE,
      color = "#63B8FF",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ pair_label,
      ncol = 1,
      scales = "free_y",
      strip.position = "left",
      labeller = label_parsed
    ) +
    scale_x_sqrt(
      limits = x_limits,
      breaks = x_breaks_use,
      labels = x_breaks_use,
      name   = "n"
    ) +
    labs(y = "", title = expression("Pooled Covariance " * s[i*j*","*"ESS"])) +
    base_theme
  
  combo <- (p_var | p_cov) + patchwork::plot_layout(widths = c(1, 1))
  
  if (!is.null(filename)) {
    ggsave(filename, combo, width = width, height = height, dpi = dpi)
  }
  
  combo
}


############################
# Skip-sampled traces (no pooling across offsets)
############################

s2skip_trace_one <- function(x, Tmax = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  N <- length(x)
  if (is.null(Tmax)) Tmax <- floor(N / 10)
  if (Tmax < 1) stop("Not enough data: N must be >= 10.")
  
  Tvec   <- seq_len(Tmax)
  np     <- N / Tvec
  s2skip <- numeric(Tmax)
  
  # T = 1
  s2skip[1] <- var(x)
  if (Tmax >= 2) {
    for (T in 2:Tmax) {
      idx <- seq.int(1L, N, by = T)
      s2skip[T] <- if (length(idx) >= 2) var(x[idx]) else NA_real_
    }
  }
  
  tibble(T = Tvec, np = np, s2skip = s2skip)
}

compute_s2skip_traces <- function(data,
                                  vars = c("co2_ppm", "humidity_pct", "soil_temp_F")) {
  stopifnot(all(vars %in% names(data)))
  
  out <- lapply(vars, function(v) {
    ans <- s2skip_trace_one(data[[v]])
    ans$variable <- v
    ans
  })
  
  bind_rows(out) |>
    mutate(variable = factor(.data$variable, levels = vars))
}

sXYskip_trace_one <- function(x, y, Tmax = NULL) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x  <- x[ok]
  y  <- y[ok]
  n  <- length(x)
  if (is.null(Tmax)) Tmax <- floor(n / 10)
  if (Tmax < 1) stop("Not enough data: n must be >= 10.")
  
  Tvec    <- seq_len(Tmax)
  np      <- n / Tvec
  sXYskip <- numeric(Tmax)
  
  # T = 1
  sXYskip[1] <- cov(x, y)
  if (Tmax >= 2) {
    for (T in 2:Tmax) {
      idx <- seq.int(1L, n, by = T)
      sXYskip[T] <- if (length(idx) >= 2) cov(x[idx], y[idx]) else NA_real_
    }
  }
  
  tibble(T = Tvec, np = np, sXYskip = sXYskip)
}

compute_pairwise_cov_skip <- function(data,
                                      vars = c("co2_ppm", "humidity_pct", "soil_temp_F")) {
  stopifnot(all(vars %in% names(data)))
  pairs <- combn(vars, 2, simplify = FALSE)
  
  out <- lapply(pairs, function(pr) {
    a <- pr[1]
    b <- pr[2]
    ans <- sXYskip_trace_one(data[[a]], data[[b]])
    ans$var_x    <- a
    ans$var_y    <- b
    ans$pair_key <- paste(a, b, sep = "_")
    ans
  })
  
  bind_rows(out)
}

############################
# Compute traces & make plots
############################

#Compute traces
vars        <- c("co2_ppm","humidity_pct","soil_temp_F")
s2p_traces  <- compute_s2p_traces(phase1, vars)
#save(s2p_traces, file = "s2p_traces.RData")
#load("s2p_traces.RData")

sxy_traces  <- compute_pairwise_cov(phase1, vars)
#save(sxy_traces, file = "sxy_traces.RData")
#load("sxy_traces.RData")

s2skip      <- compute_s2skip_traces(phase1, vars)
#save(s2skip, file = "s2skip.RData")
#load("s2skip.RData")

sXYskip     <- compute_pairwise_cov_skip(phase1, vars)
#save(sXYskip, file = "sXYskip.RData")
#load("sXYskip.RData")

# explore/choose T mask windows (T-space shading only)
plot_pooled_var_cov_T_mask_tuning(s2p_traces, sxy_traces)

plot_pooled_var_cov_mask_tuning(
  s2p_traces,
  sxy_traces,
  highlight_Ts          = c(480, 720, 960, 1440, 2160, 2880, 3600, 4320),
  highlight_T_halfwidth = c(20, 25, 20, 95, 40, 175, 60, 60)
)


# plot a single variance or covariance trace
plot_single_var_trace(s2p_traces, "co2_ppm")
plot_single_cov_trace(sxy_traces, "humidity_pct", "soil_temp_F")

# final masked figure (Figure 4)
plot_pooled_var_cov(
  s2p_traces,
  sxy_traces,
  filename = "greenhouse_var_cov.png"
)



#explore/choose mask windows
plot_pooled_var_cov_mask_tuning(s2p_traces, sxy_traces)
plot_pooled_var_cov_T_mask_tuning(s2p_traces, sxy_traces)

#plot a single variance or covariance trace
plot_single_var_trace(s2p_traces, "co2_ppm")
plot_single_cov_trace(sxy_traces, "humidity_pct", "soil_temp_F")

#plot with masking for figure 4
plot_pooled_var_cov(s2p_traces, sxy_traces,
                    filename = "greenhouse_var_cov.png")

