############################
# Choose T from pooled traces (ESS) + ACF/CCF gate
############################

# ACF / CCF gate (same as before) 
first_below_thresh_lag <- function(x, y = NULL, lag_max = 40) {
  ac <- if (is.null(y)) {
    acf(x, lag.max = lag_max, plot = FALSE, na.action = na.pass)
  } else {
    ccf(x, y, lag.max = lag_max, plot = FALSE, na.action = na.pass)
  }
  r    <- as.numeric(drop(ac$acf))
  n    <- ac$n.used
  thr  <- 2 / sqrt(n)
  lags <- as.integer(round(drop(ac$lag) * n))
  
  if (is.null(y)) {
    ok <- which(lags >= 1L & lags <= lag_max & abs(r) < thr)
    if (length(ok)) lags[ok[1]] else lag_max
  } else {
    ok <- which(abs(lags) >= 1L & abs(lags) <= lag_max & abs(r) < thr)
    if (!length(ok)) return(lag_max)
    k <- lags[ ok[ which.min(abs(lags[ok])) ] ]
    as.integer(abs(k))
  }
}

# safe rolling range (doesn't blow up when all-NA windows)
.roll_range <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  rng <- range(x, na.rm = TRUE)
  diff(rng)
}

# pick the nearest non-NA value to a target T
.value_at_T_or_nearest <- function(df, ycol, T_target) {
  y <- df[[ycol]]
  
  # If there's an exact match in T and its first y is non-NA, use that
  if (!is.na(match(T_target, df$T)) && !is.na(y[df$T == T_target][1])) {
    i <- which(df$T == T_target)[1]
    return(list(T = df$T[i], np = df$np[i], value = y[i]))
  }
  
  # Otherwise, go to the nearest T (among non-NA y's)
  ok <- which(!is.na(y))
  if (!length(ok)) {
    return(list(T = NA_integer_, np = NA_real_, value = NA_real_))
  }
  
  i <- ok[ which.min(abs(df$T[ok] - T_target)) ]
  list(T = df$T[i], np = df$np[i], value = y[i])
}

# core chooser on a single masked series
.chooseT_one <- function(df, ycol_mask,
                         method = c("acf", "acf_roll", "roll"),
                         window = 20, frac = 0.30,
                         lag_max = 40,
                         X = NULL, Y = NULL) {
  method <- match.arg(method)
  
  series <- df[[ycol_mask]]
  n      <- length(series)
  
  if (n == 0L || all(is.na(series))) {
    return(list(T_roll = NA_integer_, T_corr = NA_integer_, T_star = NA_integer_))
  }
  
  w <- min(window, n)
  
  # total range & threshold
  total_range <- diff(range(series, na.rm = TRUE))
  thr         <- frac * total_range
  
  # rolling range
  roll_rng <- rollapplyr(
    series,
    width = w,
    FUN   = .roll_range,
    partial = FALSE,
    fill    = NA_real_
  )
  
  exceed_idx <- which(roll_rng > thr)
  
  # ---- pick index for T_roll via “previous window median” rule ----
  if (length(exceed_idx)) {
    i_pick  <- exceed_idx[1]           # right edge of first window over threshold
    j_start <- i_pick - w + 1
    j_end   <- i_pick

    prev_end <- j_start - 1
    if (prev_end >= 1L) {
      prev_start <- max(1L, prev_end - w + 1L)
      idx_win    <- prev_start:prev_end
    } else {
      # no full previous window -> use earliest data up to j_end
      idx_win <- 1:j_end
    }
  } else {
    # no alarm: use the last window
    if (n <= w) {
      idx_win <- 1:n
    } else {
      idx_win <- (n - w + 1L):n
    }
  }


  window_y <- series[idx_win]
  if (all(is.na(window_y))) {
    i_roll <- idx_win[1]
  } else {
    med <- median(window_y, na.rm = TRUE)
    rel <- which(!is.na(window_y))
    k   <- rel[ which.min(abs(window_y[rel] - med)) ]
    i_roll <- idx_win[k]
  }
  
  T_roll <- df$T[i_roll]
  
  # ---- ACF/CCF gate for correlation time ----
  if (method != "roll") {
    if (is.null(Y)) {
      T_corr_floor <- first_below_thresh_lag(X, lag_max = lag_max)
    } else {
      T_corr_floor <- first_below_thresh_lag(X, y = Y, lag_max = lag_max)
    }
    T_corr <- min(max(df$T, na.rm = TRUE), T_corr_floor + 1L)
  } else {
    T_corr <- NA_integer_
  }
  
  T_star <- switch(
    method,
    "roll"     = T_roll,
    "acf"      = T_corr,
    "acf_roll" = max(T_roll, T_corr, na.rm = TRUE)
  )
  
  list(T_roll = T_roll, T_corr = T_corr, T_star = T_star)
}

# drop any rows whose T lies in one of the T0 ± halfwidth windows
drop_T_windows <- function(df,
                           T_col = "T",
                           highlight_Ts,
                           highlight_T_halfwidth) {
  if (is.null(highlight_Ts) || !length(highlight_Ts)) return(df)
  
  # recycle halfwidth if scalar
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
    if (!is.finite(T0) || !is.finite(hw)) next
    bad <- bad | (Tvals >= (T0 - hw) & Tvals <= (T0 + hw))
  }
  
  df[!bad, , drop = FALSE]
}

# MAIN: choose Ts from greenhouse traces

choose_greenhouse_Ts <- function(s2p_traces,
                                 sxy_traces,
                                 data_raw = NULL,
                                 methods = "roll",
                                 chooseT_window = 20,
                                 chooseT_frac   = 0.30,
                                 lag_max        = 40,
                                 x_limits       = c(10, 300),
                                 # T windows to exclude (usually mask_Ts)
                                 highlight_Ts          = NULL,
                                 highlight_T_halfwidth = NULL) {
  
  stopifnot(all(c("T", "np", "s2p",  "variable") %in% names(s2p_traces)))
  stopifnot(all(c("T", "np", "sXYp", "var_x", "var_y", "pair_key") %in% names(sxy_traces)))
  
  methods <- unique(methods)
  allowed <- c("acf", "acf_roll", "roll")
  if (!all(methods %in% allowed)) {
    stop("methods must be drawn from: 'acf', 'acf_roll', 'roll'.")
  }
  
  # variances: filter in n-space, then drop T-windows 
  vdat <- s2p_traces |>
    filter(np > x_limits[1], np < x_limits[2]) |>
    drop_T_windows(
      T_col               = "T",
      highlight_Ts        = highlight_Ts,
      highlight_T_halfwidth = highlight_T_halfwidth
    ) |>
    mutate(s2p_mask = s2p)
  
  # covariances: same idea 
  cdat <- sxy_traces |>
    filter(np > x_limits[1], np < x_limits[2]) |>
    drop_T_windows(
      T_col               = "T",
      highlight_Ts        = highlight_Ts,
      highlight_T_halfwidth = highlight_T_halfwidth
    ) |>
    mutate(sXYp_mask = sXYp)
  
  # choose Ts per variance 
  var_choices <- map_dfr(unique(vdat$variable), function(v) {
    dfv <- vdat |>
      filter(variable == v) |>
      arrange(T)
    
    if (!nrow(dfv)) {
      return(tibble(
        variable    = v,
        method      = methods,
        T_selected  = NA_integer_,
        np_selected = NA_real_,
        estimate    = NA_real_,
        naive_est   = NA_real_,
        T_roll      = NA_integer_,
        T_corr      = NA_integer_
      ))
    }
    
    # raw series for ACF gate
    X <- if (!is.null(data_raw)) as.numeric(data_raw[[v]]) else NULL
    if (!is.null(X)) X <- X[is.finite(X)]
    
    map_dfr(methods, function(m) {
      choice <- .chooseT_one(
        dfv, "s2p_mask",
        method = m,
        window = chooseT_window,
        frac   = chooseT_frac,
        lag_max = lag_max,
        X = X, Y = NULL
      )
      
      picked <- .value_at_T_or_nearest(dfv, "s2p_mask", choice$T_star)
      
      naive_row <- s2p_traces |>
        filter(variable == v, T == 1)
      naive_val <- if (nrow(naive_row)) naive_row$s2p[1] else NA_real_
      
      tibble(
        variable    = v,
        method      = m,
        T_selected  = as.integer(picked$T),
        np_selected = as.numeric(picked$np),
        estimate    = as.numeric(picked$value),
        naive_est   = naive_val,
        T_roll      = as.integer(choice$T_roll),
        T_corr      = as.integer(choice$T_corr)
      )
    })
  })
  
  # choose Ts per covariance pair 
  cov_choices <- map_dfr(unique(cdat$pair_key), function(pk) {
    dfp <- cdat |>
      filter(pair_key == pk) |>
      arrange(T)
    
    # If everything got dropped, still return a row per method with NAs
    if (!nrow(dfp)) {
      tmp <- sxy_traces |>
        filter(pair_key == pk) |>
        slice(1)
      vx <- tmp$var_x[1]
      vy <- tmp$var_y[1]
      return(tibble(
        var_x       = vx,
        var_y       = vy,
        pair_key    = pk,
        method      = methods,
        T_selected  = NA_integer_,
        np_selected = NA_real_,
        estimate    = NA_real_,
        naive_est   = NA_real_,
        T_roll      = NA_integer_,
        T_corr      = NA_integer_
      ))
    }
    
    vx <- dfp$var_x[1]
    vy <- dfp$var_y[1]
    
    # raw series for CCF gate
    X <- if (!is.null(data_raw)) as.numeric(data_raw[[vx]]) else NULL
    Y <- if (!is.null(data_raw)) as.numeric(data_raw[[vy]]) else NULL
    if (!is.null(X) && !is.null(Y)) {
      ok <- is.finite(X) & is.finite(Y)
      X  <- X[ok]
      Y  <- Y[ok]
    } else {
      X <- Y <- NULL
    }
    
    map_dfr(methods, function(m) {
      choice <- .chooseT_one(
        dfp, "sXYp_mask",
        method = m,
        window = chooseT_window,
        frac   = chooseT_frac,
        lag_max = lag_max,
        X = X, Y = Y
      )
      
      picked <- .value_at_T_or_nearest(dfp, "sXYp_mask", choice$T_star)
      
      naive_row <- sxy_traces |>
        filter(pair_key == pk, T == 1)
      naive_val <- if (nrow(naive_row)) naive_row$sXYp[1] else NA_real_
      
      tibble(
        var_x       = vx,
        var_y       = vy,
        pair_key    = pk,
        method      = m,
        T_selected  = as.integer(picked$T),
        np_selected = as.numeric(picked$np),
        estimate    = as.numeric(picked$value),
        naive_est   = naive_val,
        T_roll      = as.integer(choice$T_roll),
        T_corr      = as.integer(choice$T_corr)
      )
    })
  })
  
  list(
    var_choices = var_choices,
    cov_choices = cov_choices
  )
}

res_T <- choose_greenhouse_Ts(
  s2p_traces,
  sxy_traces,
  data_raw         = phase1,
  methods          = "roll",
  chooseT_window   = 30,
  chooseT_frac     = 0.2,
  x_limits         = c(10, 300),
  highlight_Ts          = mask_Ts,
  highlight_T_halfwidth = mask_T_hw
)

res_T$var_choices
res_T$cov_choices 

plot_pooled_var_cov(
  s2p_traces,
  sxy_traces,
  highlight_Ts          = mask_Ts,
  highlight_T_halfwidth = mask_T_hw,
  res_T      = res_T,
  sel_method = "roll",
  filename   = "greenhouse_var_cov_with_Tdots.png"
)
