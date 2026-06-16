rm(list=ls())
require(tidyverse)
require(lubridate)
require(ggplot2)
require(dplyr)
require(tidyr)
require(purrr)
require(patchwork)

source("colors.R")

data_dir <- "data/"

# Read all greenhouse csv files
greenhouse <- list.files(
  path = data_dir,
  pattern = "^32C4007E.*\\.csv$",  # regex pattern for file names
  full.names = TRUE
) |>
  map_dfr(~ read_csv(.x, skip = 3, col_names = c(
    "datetime_mmddyy",
    "datetime_days_since",
    "co2_ppm",
    "air_temp_F",
    "humidity_pct"
  )))

# Read all soil temperature csv files
soiltemp <- list.files(
  path = data_dir,
  pattern = "^328207FC.*\\.csv$",
  full.names = TRUE
) |>
  map_dfr(~ read_csv(.x, skip = 3, col_names = c(
    "datetime_mmddyy",
    "datetime_days_since",
    "soil_temp_F"
  )))

# Retain minute, not second
greenhouse <- greenhouse |>
  mutate(datetime_minute = floor_date(datetime_mmddyy, unit = "minute"))

soiltemp <- soiltemp |>
  mutate(datetime_minute = floor_date(datetime_mmddyy, unit = "minute"))

# Merge on the minute datetime
greenhouse <- greenhouse |>
  inner_join(
    soiltemp,
    by = "datetime_minute",
    suffix = c("_greenhouse", "_soiltemp")
  )

# Clean up 
greenhouse <- greenhouse |>
  dplyr::select(datetime_minute, co2_ppm, air_temp_F, humidity_pct, soil_temp_F) |>
  arrange(datetime_minute) |>
  rename(date=datetime_minute)

# phase 1 (not first 7 days)
phase1 <- greenhouse |>
  filter(date> ymd_hms("2025-08-12 12:00:00")) |>
  dplyr::select(-air_temp_F)

plot_greenhouse_facets <- function(data, 
                                   col = gh_cols("moody"), 
                                   red_line = NULL) {
  vars <- c("co2_ppm", "humidity_pct", "soil_temp_F")
  stopifnot(all(c("date", vars) %in% names(data)))
  if (!inherits(data$date, "POSIXct")) stop("data$date must be POSIXct")
  
  tz <- attr(data$date, "tzone")
  if (is.null(tz) || tz == ""){
    tz <- Sys.timezone()
  } 
  
  long_dat <- data |>
    select(date, all_of(vars)) |>
    pivot_longer(cols = all_of(vars), 
                 names_to = "variable", 
                 values_to = "value") |>
    mutate(variable = factor(variable, levels = vars))
  
  pretty_map <- c(
    co2_ppm      = "CO[2]~'(ppm)'",
    humidity_pct = "'Humidity (%)'",
    soil_temp_F  = "'Soil Temp ('*degree*'F)'"
  )
  
  day_starts <- seq(
    floor_date(min(long_dat$date, na.rm = TRUE), "day"),
    ceiling_date(max(long_dat$date, na.rm = TRUE), "day"),
    by = "1 day"
  )
  
  day_month_labels <- function(x) {
    x <- as.POSIXct(x, tz = tz)
    labs <- format(x, "%d", tz = tz)
    
    # Mark the first tick of each month (even if it's not the 1st because we start mid August)
    # and also mark true first-of-month ticks when present
    ym <- strftime(x, "%Y-%m", tz = tz)
    idx <- !is.na(x) & ( !duplicated(ym) | as.POSIXlt(x, tz = tz)$mday == 1 )
    
    labs[idx] <- paste0(labs[idx], "\n", strftime(x[idx], "%b", tz = tz))
    labs
  }
  
  
  p <- ggplot(long_dat, aes(x = date, y = value, group = 1)) +
    geom_line(color = col) +
    geom_vline(xintercept = day_starts, 
               linetype = "dashed", 
               color = gh_cols("grey_light")) +
    facet_wrap(
      ~ variable,
      ncol = 1,
      scales = "free_y",
      strip.position = "left",
      labeller = as_labeller(pretty_map, default = label_parsed)
    ) +
    scale_x_datetime(
      breaks = day_starts,
      labels = day_month_labels,
      expand = expansion(mult = c(0, 0))
    ) +
    guides(x = guide_axis(check.overlap = TRUE)) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "plain"),
      axis.text.x = element_text(margin = margin(t = 6))
    ) +
    coord_cartesian(clip = "on") +
    theme(plot.margin = margin(b = 12, r = 6, l = 6, t = 6))
  
  if (!is.null(red_line)) {
    p <- p + geom_vline(xintercept = as.POSIXct(red_line, tz = tz), color = "red")
  }
  
  p
}


greenhousetrace <- plot_greenhouse_facets(phase1, col = gh_cols("moody"))
ggsave("greenhousetrace.pdf", greenhousetrace)

##############
#make acfs and ccfs
require(dplyr)
require(tibble)
require(purrr)
require(ggplot2)
require(patchwork)  

# settings
lag_max <- 20000
period  <- 1440  # 24 hours in minutes

# Labels (row order)
acf_labels <- c(
  co2_ppm      = "CO₂",
  humidity_pct = "Humidity",
  soil_temp_F  = "Soil Temp"
)

ccf_pairs <- list(
  c("co2_ppm",      "humidity_pct"),
  c("humidity_pct", "soil_temp_F"),
  c("soil_temp_F",  "co2_ppm")
)
ccf_labels <- c(
  "co2_ppm|humidity_pct"      = "CO₂ vs Humidity",
  "humidity_pct|soil_temp_F"  = "Humidity vs Soil Temp",
  "soil_temp_F|co2_ppm"       = "Soil Temp vs CO₂"
)

# labels as plotmath strings
acf_labels <- c(
  co2_ppm      = "CO[2]",
  humidity_pct = "Humidity",
  soil_temp_F  = "Soil~Temp"
)

ccf_labels <- c(
  "co2_ppm|humidity_pct"      = "CO[2]~'vs'~Humidity",
  "humidity_pct|soil_temp_F"  = "Humidity~'vs'~'Soil Temp'",
  "soil_temp_F|co2_ppm"       = "'Soil Temp'~'vs'~CO[2]"
)


# Maintain consistent row order
acf_levels <- unname(acf_labels)
ccf_levels <- ccf_labels[c("co2_ppm|humidity_pct",
                           "humidity_pct|soil_temp_F",
                           "soil_temp_F|co2_ppm")] |>
  unname()

# Helpers
acf_df <- function(x, lag_max) {
  a <- stats::acf(x, lag.max = lag_max, plot = FALSE, na.action = na.pass)
  tibble(
    lag  = as.integer(a$lag),
    corr = as.numeric(a$acf)
  ) |> dplyr::filter(lag >= 0)     # ACF: nonnegative lags only
}

ccf_df <- function(x, y, lag_max) {
  c <- stats::ccf(x, y, lag.max = lag_max, plot = FALSE, na.action = na.pass)
  tibble(
    lag  = as.integer(c$lag),      # CCF: keep ± lags
    corr = as.numeric(c$acf)
  )
}

# Build tidy ACF for the left column
acfs <- purrr::imap_dfr(names(acf_labels), function(nm, i) {
  out <- acf_df(phase1[[nm]], lag_max)
  out$panel_label <- acf_labels[[nm]]
  out
}) |>
  dplyr::mutate(panel_label = factor(panel_label, levels = acf_levels))

vlines_acf <- tibble(lag = seq(0, lag_max, by = period))

p_acf <- ggplot(acfs, aes(lag, corr)) +
  geom_segment(aes(xend = lag, y = 0, yend = corr), 
               linewidth = 0.35 , 
               col = gh_cols("moody")) +
  geom_hline(yintercept = 0, 
             linewidth = 0.3,
             color = gh_cols("grey_dark")) +
  geom_vline(data = vlines_acf, 
             aes(xintercept = lag),
             inherit.aes = FALSE, 
             linewidth = 0.3, 
             color = gh_cols("grey_mid"),
             linetype = "longdash", 
             alpha = 0.5) +
  coord_cartesian(xlim = c(0, lag_max), 
                  ylim = c(-1, 1)) +   
  facet_wrap(~ panel_label, 
             ncol = 1, 
             scales = "fixed", 
             strip.position = "top") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Lag", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position      = "none",
    strip.text           = element_text(face = "bold"),
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank()    # no vertical grid lines
  )

# Build tidy CCF for right column
ccfs <- purrr::imap_dfr(ccf_pairs, function(v, i) {
  xnm <- v[1]; ynm <- v[2]
  key <- paste0(xnm, "|", ynm)
  out <- ccf_df(phase1[[xnm]], phase1[[ynm]], lag_max)
  out$panel_label <- ccf_labels[[key]]
  out
}) |>
  dplyr::mutate(panel_label = factor(panel_label, levels = ccf_levels))

vlines_ccf <- tibble(lag = seq(-lag_max, lag_max, by = period))

p_ccf <- ggplot(ccfs, aes(lag, corr)) +
  geom_segment(aes(xend = lag, y = 0, yend = corr), 
               linewidth = 0.35, 
               col = gh_cols("moody")) +
  geom_hline(yintercept = 0, 
             linewidth = 0.3, 
             color = gh_cols("grey_dark")) +
  geom_vline(data = vlines_ccf, 
             aes(xintercept = lag),
             inherit.aes = FALSE, 
             linewidth = 0.3, 
             color = gh_cols("grey_mid"),
             linetype = "longdash", 
             alpha = 0.5) +
  coord_cartesian(xlim = c(-lag_max, lag_max), 
                  ylim = c(-1, 1)) +
  facet_wrap(~ panel_label, 
             ncol = 1,
             scales = "fixed", 
             strip.position = "top") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Lag", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position      = "none",
    strip.text           = element_text(face = "bold"),
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank()
  )

# ---------------------------
# Combine: 3 rows x 2 columns
# ---------------------------
# when you facet, add labeller = label_parsed
p_acf <- p_acf + facet_wrap(~ panel_label, 
                            ncol = 1, 
                            strip.position = "top",
                            labeller = label_parsed)
p_ccf <- p_ccf + facet_wrap(~ panel_label, 
                            ncol = 1, 
                            strip.position = "top",
                            labeller = label_parsed)

greenhouseacfccfs <- p_acf | p_ccf

ggsave("acfccf.pdf", greenhouseacfccfs)

############################
# Plot the pooled traces

# Minimal deps
library(dplyr)
library(ggplot2)
library(patchwork)

# ---------- pooled variance traces (with progress) ----------
s2p_trace_one <- function(x, Tmax = NULL, progress = TRUE) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  N <- length(x)
  if (is.null(Tmax)) Tmax <- floor(N/10)
  if (Tmax < 1) stop("Not enough data: N must be >= 10.")
  
  Tvec <- seq_len(Tmax)
  np   <- N / Tvec                      # <-- continuous, not floor()
  s2p  <- numeric(Tmax)
  
  if (progress) {
    message(sprintf("[variance] N=%d, Tmax=%d", N, Tmax)); flush.console()
    pb <- utils::txtProgressBar(min = 1, max = Tmax, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  
  # T = 1
  s2p[1] <- stats::var(x); if (progress) utils::setTxtProgressBar(pb, 1L)
  
  if (Tmax >= 2) {
    for (T in 2:Tmax) {
      s2_t <- numeric(T)
      for (i in 1:T) {
        idx <- seq.int(i, N, by = T)
        s2_t[i] <- if (length(idx) >= 2) stats::var(x[idx]) else NA_real_
      }
      s2p[T] <- mean(s2_t, na.rm = TRUE)
      if (progress) utils::setTxtProgressBar(pb, T)
    }
  }
  tibble::tibble(T = Tvec, np = np, s2p = s2p)
}


compute_s2p_traces <- function(data,
                               vars = c("co2_ppm","humidity_pct","soil_temp_F"),
                               progress = TRUE) {
  stopifnot(all(vars %in% names(data)))
  if (progress) {
    message(sprintf("Computing pooled variances for %d variables ...", length(vars)))
    flush.console()
  }
  out <- lapply(seq_along(vars), function(j) {
    v <- vars[j]
    if (progress) { message(sprintf("  [%d/%d] %s", j, length(vars), v)); flush.console() }
    ans <- s2p_trace_one(data[[v]], progress = progress)
    ans$variable <- v
    ans
  })
  dplyr::bind_rows(out) |>
    dplyr::mutate(variable = factor(.data$variable, levels = vars))
}

# ---------- pooled covariance traces (with progress) ----------
compute_pairwise_cov <- function(data,
                                 vars = c("co2_ppm","humidity_pct","soil_temp_F"),
                                 progress = TRUE) {
  stopifnot(all(vars %in% names(data)))
  pairs <- combn(vars, 2, simplify = FALSE)
  
  lapply(seq_along(pairs), function(k) {
    pr <- pairs[[k]]; a <- pr[1]; b <- pr[2]
    if (progress) message(sprintf("  [%d/%d] %s × %s", k, length(pairs), a, b))
    
    x <- as.numeric(data[[a]]); y <- as.numeric(data[[b]])
    ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
    n <- length(x); if (n < 2) return(NULL)
    
    Tmax <- floor(n / 10)
    Tvec <- 1:Tmax
    np   <- n / Tvec                 # <-- use n/T (continuous)
    sXYp <- numeric(Tmax)
    
    if (progress) {
      pb <- utils::txtProgressBar(min = 1, max = Tmax, style = 3)
      on.exit(close(pb), add = TRUE)
    }
    
    # T = 1
    sXYp[1] <- stats::cov(x, y)
    if (progress) utils::setTxtProgressBar(pb, 1L)
    
    if (Tmax >= 2) {
      for (T in 2:Tmax) {
        cov_t <- numeric(T)
        for (i in 1:T) {
          idx <- seq.int(i, n, by = T)
          cov_t[i] <- if (length(idx) >= 2) stats::cov(x[idx], y[idx]) else NA_real_
        }
        sXYp[T] <- mean(cov_t, na.rm = TRUE)
        if (progress) utils::setTxtProgressBar(pb, T)
      }
    }
    
    tibble::tibble(
      var_x = a, var_y = b, pair_key = paste(a,b,sep="_"),
      T = Tvec, np = np, sXYp = sXYp
    )
  }) |>
    dplyr::bind_rows()
}

# ---------- plotting: left = variances, right = covariances ----------
# helper from earlier (unchanged)

# helpers
pretty_var <- function(v){
  c(co2_ppm="CO[2]", humidity_pct="Humidity", soil_temp_F="Soil~Temperature")[v]
}
# stacked, two-line facet label for covariance: atop(top, bottom)
pair_label_stacked <- function(x, y){
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
  bad <- rep(FALSE, length(Tvals))
  
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

plot_pooled_var_cov <- function(
    var_traces,
    cov_traces,
    col = if (exists("gh_cols")) gh_cols("moody") else "#36648B",
    filename = NULL,
    width = 10, height = 8, dpi = 300,
    x_limits = c(10, 300),
    x_breaks  = c(10, 50, 100, 150, 200, 250, 300),
    
    # --- the ONLY masking: these T windows ---
    highlight_Ts          = c(480, 720, 960, 1440, 2160, 2880, 3600, 4320),
    highlight_T_halfwidth = c(20,  25,  20,   95,   40,   175,  60,   60),
    
    # --- show selected T's as diamonds ---
    res_T          = NULL,          # list(var_choices, cov_choices) from choose_greenhouse_Ts
    sel_method     = "roll",        # which method's T to show ("roll", "acf", "acf_roll")
    sel_point_col  = if (exists("gh_cols")) gh_cols("orange") else "#228B22"
) {
  stopifnot(all(c("T","np","s2p","variable") %in% names(var_traces)))
  stopifnot(all(c("T","np","sXYp","var_x","var_y","pair_key") %in% names(cov_traces)))
  
  # Filter to the n-range you care about, but keep full T info
  vdat <- var_traces |>
    dplyr::filter(np > x_limits[1], np < x_limits[2]) |>
    mask_T_windows(
      T_col = "T",
      y_col = "s2p",
      highlight_Ts          = highlight_Ts,
      highlight_T_halfwidth = highlight_T_halfwidth
    )
  
  cdat <- cov_traces |>
    dplyr::mutate(pair_label = pair_label_stacked(var_x, var_y)) |>
    dplyr::filter(np > x_limits[1], np < x_limits[2]) |>
    mask_T_windows(
      T_col = "T",
      y_col = "sXYp",
      highlight_Ts          = highlight_Ts,
      highlight_T_halfwidth = highlight_T_halfwidth
    )
  
  # T = 1 references for the dashed lines (use *unmasked* traces)
  t1_var <- var_traces |>
    dplyr::filter(T == 1) |>
    dplyr::select(variable, s2p)
  
  t1_cov <- cov_traces |>
    dplyr::filter(T == 1) |>
    dplyr::mutate(pair_label = pair_label_stacked(var_x, var_y)) |>
    dplyr::select(pair_label, sXYp)
  
  # ------- build data for selected T points from res_T (if supplied) -------
  var_pts <- NULL
  cov_pts <- NULL
  
  if (!is.null(res_T)) {
    # variance selections
    if ("var_choices" %in% names(res_T)) {
      var_pts <- res_T$var_choices |>
        dplyr::filter(method == sel_method,
                      !is.na(np_selected),
                      !is.na(estimate)) |>
        dplyr::select(variable, np_selected, estimate) |>
        dplyr::distinct() |>
        dplyr::filter(np_selected > x_limits[1], np_selected < x_limits[2])
      
      if (!is.null(vdat$variable)) {
        var_pts$variable <- factor(var_pts$variable,
                                   levels = levels(vdat$variable))
      }
    }
    
    # covariance selections
    if ("cov_choices" %in% names(res_T)) {
      cov_pts <- res_T$cov_choices |>
        dplyr::filter(method == sel_method,
                      !is.na(np_selected),
                      !is.na(estimate)) |>
        dplyr::mutate(pair_label = pair_label_stacked(var_x, var_y)) |>
        dplyr::select(pair_label, np_selected, estimate) |>
        dplyr::distinct() |>
        dplyr::filter(np_selected > x_limits[1], np_selected < x_limits[2])
    }
  }
  
  # Labels as in your original function
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
  
  # ---- left: pooled variance (masked only in the specified T ranges) ----
  p_var <- ggplot(vdat, aes(x = np, y = s2p_mask)) +
    geom_line(color = col, na.rm = TRUE) +
    # green diamonds at selected T/n,estimate
    { if (!is.null(var_pts) && nrow(var_pts) > 0)
      geom_point(
        data = var_pts,
        inherit.aes = FALSE,
        aes(x = np_selected, y = estimate),
        shape = 23,              # filled diamond
        size  = 2.8,
        stroke = 0.4,
        color = "black",
        fill  = sel_point_col
      )
    } +
    geom_hline(
      data = t1_var, aes(yintercept = s2p),
      inherit.aes = FALSE,
      color = "#63B8FF",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ variable, ncol = 1, scales = "free_y",
      strip.position = "left",
      labeller = var_labeller
    ) +
    scale_x_sqrt(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_breaks,
      name   = expression(n[t])
    ) +
    labs(y = "", title = expression("Pooled Variance " * s[i*","*"ESS"]^2)) +
    base_theme
  
  # ---- right: pooled covariance (masked only in the specified T ranges) ----
  p_cov <- ggplot(cdat, aes(x = np, y = sXYp_mask)) +
    geom_line(color = col, na.rm = TRUE) +
    { if (!is.null(cov_pts) && nrow(cov_pts) > 0)
      geom_point(
        data = cov_pts,
        inherit.aes = FALSE,
        aes(x = np_selected, y = estimate),
        shape = 23,
        size  = 2.8,
        stroke = 0.4,
        color = "black",
        fill  = sel_point_col
      )
    } +
    geom_hline(
      data = t1_cov, aes(yintercept = sXYp),
      inherit.aes = FALSE,
      color = "#63B8FF",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ pair_label, ncol = 1, scales = "free_y",
      strip.position = "left",
      labeller = label_parsed
    ) +
    scale_x_sqrt(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_breaks,
      name   = expression(n[t])
    ) +
    labs(y = "", title = expression("Pooled Covariance " * s[i*j*","*"ESS"])) +
    base_theme
  
  combo <- (p_var | p_cov) + patchwork::plot_layout(widths = c(1, 1))
  
  if (!is.null(filename)) {
    ggsave(filename, combo, width = width, height = height, dpi = dpi)
  }
  
  combo
}



# ---------- run (example with your phase1; N = nrow(phase1) = 41577) ----------
vars <- c("co2_ppm","humidity_pct","soil_temp_F")

#s2p_traces  <- compute_s2p_traces(phase1, vars)         # variances, T=1..floor(N/10)
#save(s2p_traces, file = "s2p_traces_full.RData")
#load("s2p_traces_full.RData")
#sxy_traces  <- compute_pairwise_cov(phase1, vars)       # covariances, T=1..floor(N/10)
#save(sxy_traces, file = "sxy_traces_full.RData")
#load("sxy_traces_full.RData")
# build the usual pooled traces (you already have these)

plot_pooled_var_cov(
  s2p_traces,
  sxy_traces,
  filename = "greenhouse_var_cov.png"
)


###############
# masking debug
##############
plot_pooled_var_cov_debug <- function(
    var_traces,
    cov_traces,
    col = if (exists("gh_cols")) gh_cols("moody") else "#36648B",
    filename = NULL,
    width = 10, height = 8, dpi = 300,
    x_limits = c(10, 300),
    x_breaks  = c(10, 50, 100, 150, 200, 250, 300),
    
    # --- shading controls ---
    # 1) resonance-like periods (like your current mask_resonances)
    highlight_periods    = c(720, 1440, 2160, 2880),  # 12h, 24h, 36h, 48h
    highlight_np_halfwidths = c(0.8, 1.0, 0.3, 0.45),
    use_period_boxes     = TRUE,
    
    # 2) explicit skip T windows (e.g. T = 60 ± 5)
    #    list(list(T = 60, halfwidth = 5), list(T = 120, halfwidth = 10))
    highlight_T_windows  = NULL,
    
    # 3) explicit n windows (e.g. n = 11.5 ± 0.5)
    #    list(list(n = 11.5, halfwidth = 0.5), list(n = 31, halfwidth = 1))
    highlight_n_windows  = NULL,
    
    shade_color = "hotpink",
    shade_alpha = 0.12,
    
    # more ticks for exploring n
    dense_ticks = FALSE,
    n_dense_breaks = 12
) {
  stopifnot(all(c("T","np","s2p","variable") %in% names(var_traces)))
  stopifnot(all(c("T","np","sXYp","var_x","var_y","pair_key") %in% names(cov_traces)))
  
  # keep original curves, just limit in n
  vdat <- var_traces |>
    dplyr::filter(np > x_limits[1], np < x_limits[2])
  
  cdat <- cov_traces |>
    dplyr::filter(np > x_limits[1], np < x_limits[2]) |>
    dplyr::mutate(pair_label = pair_label_stacked(var_x, var_y))
  
  # estimate N from traces (as you do elsewhere: np * T ≈ N)
  N_est <- round(stats::median(var_traces$np * var_traces$T, na.rm = TRUE))
  
  # -------------------------
  # build highlight rectangles in n-space
  # -------------------------
  rects <- tibble::tibble(
    xmin   = numeric(0),
    xmax   = numeric(0),
    source = character(0),
    label  = character(0)
  )
  
  # 1) resonance-like periods → centers in np
  if (use_period_boxes && !is.null(highlight_periods) && length(highlight_periods) > 0) {
    centers <- N_est / highlight_periods
    # recycle halfwidths if scalar
    if (length(highlight_np_halfwidths) == 1L &&
        length(centers) > 1L) {
      highlight_np_halfwidths <- rep(highlight_np_halfwidths, length(centers))
    }
    stopifnot(length(centers) == length(highlight_np_halfwidths))
    
    rects_period <- tibble::tibble(
      xmin   = centers - highlight_np_halfwidths,
      xmax   = centers + highlight_np_halfwidths,
      source = "period",
      label  = paste0("period=", highlight_periods)
    )
    rects <- dplyr::bind_rows(rects, rects_period)
  }
  
  # 2) skip T windows → convert to n ranges: n = N / T
  if (!is.null(highlight_T_windows) && length(highlight_T_windows) > 0) {
    rects_T <- purrr::map_dfr(highlight_T_windows, function(w) {
      T0 <- w$T
      hw <- w$halfwidth
      T_low  <- pmax(1, T0 - hw)
      T_high <- pmax(T_low + 1e-8, T0 + hw)
      n_low  <- N_est / T_high
      n_high <- N_est / T_low
      tibble::tibble(
        xmin   = pmin(n_low, n_high),
        xmax   = pmax(n_low, n_high),
        source = "T",
        label  = paste0("T≈", T0)
      )
    })
    rects <- dplyr::bind_rows(rects, rects_T)
  }
  
  # 3) direct n windows
  if (!is.null(highlight_n_windows) && length(highlight_n_windows) > 0) {
    rects_n <- purrr::map_dfr(highlight_n_windows, function(w) {
      n0 <- w$n
      hw <- w$halfwidth
      tibble::tibble(
        xmin   = n0 - hw,
        xmax   = n0 + hw,
        source = "n",
        label  = paste0("n≈", n0)
      )
    })
    rects <- dplyr::bind_rows(rects, rects_n)
  }
  
  # clip to plotting window & drop degenerate ranges
  if (nrow(rects) > 0) {
    rects <- rects |>
      dplyr::mutate(
        xmin = pmax(xmin, x_limits[1]),
        xmax = pmin(xmax, x_limits[2])
      ) |>
      dplyr::filter(xmin < xmax)
  }
  
  # -------------------------
  # T = 1 reference lines
  # -------------------------
  t1_var <- var_traces |>
    dplyr::filter(T == 1) |>
    dplyr::select(variable, s2p)
  
  t1_cov <- cov_traces |>
    dplyr::filter(T == 1) |>
    dplyr::mutate(pair_label = pair_label_stacked(var_x, var_y)) |>
    dplyr::select(pair_label, sXYp)
  
  # tick strategy
  x_breaks_use <- if (dense_ticks) {
    pretty(x_limits, n = n_dense_breaks)
  } else {
    x_breaks
  }
  
  # faceting labels (same as your plot_pooled_var_cov)
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
  
  # -------------------------
  # left: pooled variance (full curve)
  # -------------------------
  p_var <- ggplot(vdat, aes(x = np, y = s2p)) +
    { if (nrow(rects) > 0)
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
      data = t1_var, aes(yintercept = s2p),
      inherit.aes = FALSE, color = "#63B8FF",
      linewidth = 0.4, linetype = "dashed"
    ) +
    facet_wrap(
      ~ variable, ncol = 1, scales = "free_y",
      strip.position = "left",
      labeller = var_labeller
    ) +
    scale_x_sqrt(
      limits = x_limits,
      breaks = x_breaks_use,
      labels = x_breaks_use,
      name   = "n"
    ) +
    labs(y = "", title = expression("Pooled Variance " * s[i*","*"ESS"]^2)) +
    base_theme
  
  # -------------------------
  # right: pooled covariance (full curve)
  # -------------------------
  p_cov <- ggplot(cdat, aes(x = np, y = sXYp)) +
    { if (nrow(rects) > 0)
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
      data = t1_cov, aes(yintercept = sXYp),
      inherit.aes = FALSE, color = "#63B8FF",
      linewidth = 0.4, linetype = "dashed"
    ) +
    facet_wrap(
      ~ pair_label, ncol = 1, scales = "free_y",
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

plot_pooled_var_cov_debug(
  s2p_traces,
  sxy_traces,
  highlight_periods      = c(720, 1440, 2160, 2880, 3600, 4320),
  highlight_np_halfwidths = c(0.8, 1.0, 0.3, 0.45, .4, .4),
  use_period_boxes       = TRUE,
  dense_ticks            = TRUE, 
  filename = "greenhouse_var_cov_debug.png"
)

plot_pooled_var_cov_T_debug <- function(
    var_traces,
    cov_traces,
    col = if (exists("gh_cols")) gh_cols("moody") else "#36648B",
    filename = NULL,
    width = 10, height = 8, dpi = 300,
    
    # --- x-axis in T ---
    T_limits = c(10, 5000),
    T_breaks = c(10, 100, 500, 1000, 2000, 4000),
    log10_x  = TRUE,      # set FALSE if you want linear T
    
    # --- resonance T's (centers) and windows around them ---
    # centers in T (minutes), e.g. 12h, 24h, 36h, 48h, 60h, 72h
    highlight_Ts          = c(720, 1440, 2160, 2880, 3600, 4320),
    highlight_T_halfwidth = 100,    # can be scalar OR same length as highlight_Ts
    
    shade_color = "hotpink",
    shade_alpha = 0.15
) {
  stopifnot(all(c("T","np","s2p","variable") %in% names(var_traces)))
  stopifnot(all(c("T","np","sXYp","var_x","var_y","pair_key") %in% names(cov_traces)))
  
  # --- keep full curves, just restrict T-range ---
  vdat <- var_traces |>
    dplyr::filter(T >= T_limits[1], T <= T_limits[2])
  
  cdat <- cov_traces |>
    dplyr::filter(T >= T_limits[1], T <= T_limits[2]) |>
    dplyr::mutate(pair_label = pair_label_stacked(var_x, var_y))
  
  # --- T = 1 references (same idea as before, but now w.r.t T) ---
  t1_var <- var_traces |>
    dplyr::filter(T == 1) |>
    dplyr::select(variable, s2p)
  
  t1_cov <- cov_traces |>
    dplyr::filter(T == 1) |>
    dplyr::mutate(pair_label = pair_label_stacked(var_x, var_y)) |>
    dplyr::select(pair_label, sXYp)
  
  # --- build T-windows as rectangles ---
  if (length(highlight_Ts) > 0) {
    # recycle halfwidths if scalar
    if (length(highlight_T_halfwidth) == 1L &&
        length(highlight_Ts) > 1L) {
      highlight_T_halfwidth <- rep(highlight_T_halfwidth, length(highlight_Ts))
    }
    stopifnot(length(highlight_Ts) == length(highlight_T_halfwidth))
    
    rects <- tibble::tibble(
      Tmin  = highlight_Ts - highlight_T_halfwidth,
      Tmax  = highlight_Ts + highlight_T_halfwidth,
      label = paste0("T≈", highlight_Ts)
    ) |>
      dplyr::mutate(
        Tmin = pmax(Tmin, T_limits[1]),
        Tmax = pmin(Tmax, T_limits[2])
      ) |>
      dplyr::filter(Tmin < Tmax)
  } else {
    rects <- tibble::tibble(Tmin = numeric(0), Tmax = numeric(0), label = character(0))
  }
  
  # --- labels like your existing plot_pooled_var_cov ---
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
  
  # choose x scale
  x_scale_fun <- if (log10_x) {
    scale_x_log10(
      limits = T_limits,
      breaks = T_breaks,
      labels = T_breaks,
      name   = "T (skip size, minutes)"
    )
  } else {
    scale_x_continuous(
      limits = T_limits,
      breaks = T_breaks,
      labels = T_breaks,
      name   = "T (skip size, minutes)"
    )
  }
  
  # --- left: pooled variance vs T, full curve + shaded windows ---
  p_var <- ggplot(vdat, aes(x = T, y = s2p)) +
    { if (nrow(rects) > 0)
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
      data = t1_var, aes(yintercept = s2p),
      inherit.aes = FALSE, color = "#63B8FF",
      linewidth = 0.4, linetype = "dashed"
    ) +
    facet_wrap(
      ~ variable, ncol = 1, scales = "free_y",
      strip.position = "left",
      labeller = var_labeller
    ) +
    x_scale_fun +
    labs(y = "", title = expression("Pooled Variance " * s[i*","*"ESS"]^2)) +
    base_theme
  
  # --- right: pooled covariance vs T, full curve + shaded windows ---
  p_cov <- ggplot(cdat, aes(x = T, y = sXYp)) +
    { if (nrow(rects) > 0)
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
      data = t1_cov, aes(yintercept = sXYp),
      inherit.aes = FALSE, color = "#63B8FF",
      linewidth = 0.4, linetype = "dashed"
    ) +
    facet_wrap(
      ~ pair_label, ncol = 1, scales = "free_y",
      strip.position = "left",
      labeller = label_parsed
    ) +
    x_scale_fun +
    labs(y = "", title = expression("Pooled Covariance " * s[i*j*","*"ESS"])) +
    base_theme
  
  combo <- (p_var | p_cov) + patchwork::plot_layout(widths = c(1, 1))
  
  if (!is.null(filename)) {
    ggsave(filename, combo, width = width, height = height, dpi = dpi)
  }
  
  combo
}


plot_pooled_var_cov_T_debug(
  s2p_traces,
  sxy_traces,
  highlight_Ts          = c(480, 720, 960, 1440, 2160, 2880, 3600, 4320),
  #8, 12, 16, 24, 36, 48, 60, 72 hours
  highlight_T_halfwidth = c(20, 25, 20, 65,  40,   110,   60,   60),
  log10_x               = FALSE,
  filename="greenhouse_var_cov_debug.png"
)


# Now look at the ESS versus skip
# do the same for skip samples: 
# ---------- skip-sampled traces (no pooling across offsets) ----------
s2skip_trace_one <- function(x, Tmax = NULL) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  N <- length(x)
  if (is.null(Tmax)) Tmax <- floor(N/10)
  if (Tmax < 1) stop("Not enough data: N must be >= 10.")
  Tvec   <- seq_len(Tmax)
  np     <- N / Tvec
  s2skip <- numeric(Tmax)
  # T = 1
  s2skip[1] <- stats::var(x)
  if (Tmax >= 2) {
    for (T in 2:Tmax) {
      idx <- seq.int(1L, N, by = T)
      s2skip[T] <- if (length(idx) >= 2) stats::var(x[idx]) else NA_real_
    }
  }
  tibble::tibble(T = Tvec, np = np, s2skip = s2skip)
}

compute_s2skip_traces <- function(data,
                                  vars = c("co2_ppm","humidity_pct","soil_temp_F")) {
  stopifnot(all(vars %in% names(data)))
  out <- lapply(vars, function(v) {
    ans <- s2skip_trace_one(data[[v]])
    ans$variable <- v
    ans
  })
  dplyr::bind_rows(out) |>
    dplyr::mutate(variable = factor(.data$variable, levels = vars))
}

sXYskip_trace_one <- function(x, y, Tmax = NULL) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (is.null(Tmax)) Tmax <- floor(n/10)
  if (Tmax < 1) stop("Not enough data: n must be >= 10.")
  Tvec    <- seq_len(Tmax)
  np      <- n / Tvec
  sXYskip <- numeric(Tmax)
  # T = 1
  sXYskip[1] <- stats::cov(x, y)
  if (Tmax >= 2) {
    for (T in 2:Tmax) {
      idx <- seq.int(1L, n, by = T)
      sXYskip[T] <- if (length(idx) >= 2) stats::cov(x[idx], y[idx]) else NA_real_
    }
  }
  tibble::tibble(T = Tvec, np = np, sXYskip = sXYskip)
}

compute_pairwise_cov_skip <- function(data,
                                      vars = c("co2_ppm","humidity_pct","soil_temp_F")) {
  stopifnot(all(vars %in% names(data)))
  pairs <- combn(vars, 2, simplify = FALSE)
  dplyr::bind_rows(lapply(pairs, function(pr) {
    a <- pr[1]; b <- pr[2]
    ans <- sXYskip_trace_one(data[[a]], data[[b]])
    ans$var_x <- a; ans$var_y <- b; ans$pair_key <- paste(a,b,sep="_")
    ans
  }))
}

# ========== Single-series plotter: ESS vs Skip vs Naive (variance or covariance) ==========
# Needs:
#   s2p_traces, sxy_traces (ESS/pooled traces)
#   s2skip_traces, sXYskip_traces (skip traces)
#   data_raw: the raw data frame with columns for variables to compute naive from
#
# Styling:
#   - Skip: solid light grey
#   - ESS : solid accent
#   - Naive: dashed horizontal line, included in legend

plot_single_ess_skip_naive <- function(
    select,                      # "co2_ppm"  OR  c("co2_ppm","soil_temp_F")
    data_raw,                    # data.frame with those columns; used for naive
    s2p_traces = NULL, sxy_traces = NULL,
    s2skip_traces = NULL, sXYskip_traces = NULL,
    type = c("auto","var","cov"),
    # style
    col_ess  = if (exists("gh_cols")) gh_cols("moody") else "#36648B",
    col_skip = if (exists("gh_cols")) gh_cols("grey_light") else "#C8C8C8",
    col_naive = "#63B8FF",
    lwd_ess  = .7, lwd_skip = 0.7, lwd_naive = 0.5,
    # x-axis
    x_limits = c(10, 300),
    x_breaks = c(10, 50, 100, 150, 200, 250, 300),
    # masking (old NA-masking, still available if you want it)
    mask_pooled = FALSE,
    mask_skip   = FALSE,
    # manual windows only used if masking is on (old n-based stuff)
    cov_windows = list(
      "co2_ppm|humidity_pct"      = list(n = c(42, 70, 85, 97), tol = c(0.5, 3, 2, 1)),
      "soil_temp_F|co2_ppm"       = list(n = c(42, 70, 95),     tol = c(0.5, 5, 0.3)),
      "humidity_pct|soil_temp_F"  = list(n = c(11.5, 32),       tol = c(0.5, 1.0))
    ),
    var_manual = list(           # only used if mask_* = TRUE
      hum_n  = 11.5,  hum_tol  = 0.5,
      soil_n = c(11.5, 32), soil_tol = c(0.5, 1.0),
      co2_n  = numeric(0),  co2_tol  = numeric(0)
    ),
    # --- NEW: T-based shaded regions + labels (visual only) ---
    show_T_boxes          = FALSE,
    highlight_Ts          = c(480, 720, 960, 1440, 2160, 2880),
    highlight_T_halfwidth = c(8.2,   15,  23.1,   65,   78,   183.6),
    highlight_labels      = NULL,     # if NULL, we use hours like "8h", "12h", ...
    shade_alpha           = 0.15,
    label_T_boxes         = TRUE,
    label_position        = c("top","bottom"),
    # output
    title_prefix = NULL,
    filename = NULL, width = 7, height = 4, dpi = 300
) {
  stopifnot(is.data.frame(data_raw))
  type <- match.arg(type)
  label_position <- match.arg(label_position)
  if (type == "auto") type <- if (length(select) == 1) "var" else "cov"
  
  # small label helpers if yours aren't in scope
  .pretty_var <- if (exists("pretty_var")) pretty_var else function(v){
    c(co2_ppm="CO[2]", humidity_pct="Humidity", soil_temp_F="Soil~Temp")[v]
  }
  .pair_label <- if (exists("pair_label_stacked")) pair_label_stacked else function(x, y){
    paste0("atop(", .pretty_var(x), ", ", .pretty_var(y), ")")
  }
  
  # helper: build (np) rectangles & labels from T windows using ESS grid
  build_T_rects_np <- function(df,
                               T_col = "T", np_col = "np",
                               highlight_Ts, highlight_T_halfwidth,
                               highlight_labels, x_limits) {
    if (!show_T_boxes || is.null(highlight_Ts) || length(highlight_Ts) == 0) {
      return(NULL)
    }
    # recycle halfwidth if scalar
    if (length(highlight_T_halfwidth) == 1L && length(highlight_Ts) > 1L) {
      highlight_T_halfwidth <- rep(highlight_T_halfwidth, length(highlight_Ts))
    }
    stopifnot(length(highlight_Ts) == length(highlight_T_halfwidth))
    
    # default labels: convert T (minutes) to "xh"
    if (is.null(highlight_labels)) {
      hrs <- highlight_Ts / 60
      # if integer-ish hours, print as integer; else one decimal
      hrs_round <- ifelse(abs(hrs - round(hrs)) < 1e-6,
                          as.character(round(hrs)),
                          format(round(hrs, 1), trim = TRUE))
      highlight_labels <- paste0(hrs_round, "h")
    } else {
      stopifnot(length(highlight_labels) == length(highlight_Ts))
    }
    
    rects <- purrr::map_dfr(seq_along(highlight_Ts), function(k) {
      T0 <- highlight_Ts[k]
      hw <- highlight_T_halfwidth[k]
      lo <- T0 - hw
      hi <- T0 + hw
      sub <- df[df[[T_col]] >= lo & df[[T_col]] <= hi, , drop = FALSE]
      if (!nrow(sub)) return(NULL)
      x_min <- max(min(sub[[np_col]], na.rm = TRUE), x_limits[1])
      x_max <- min(max(sub[[np_col]], na.rm = TRUE), x_limits[2])
      if (!is.finite(x_min) || !is.finite(x_max) || x_min >= x_max) return(NULL)
      tibble::tibble(
        xmin   = x_min,
        xmax   = x_max,
        x_label = (x_min + x_max) / 2,
        label  = highlight_labels[k]
      )
    })
    if (!nrow(rects)) return(NULL)
    dplyr::distinct(rects)
  }
  
  # position for text labels
  y_pos    <- if (label_position == "top")  Inf else -Inf
  vjust_lab <- if (label_position == "top") 1.2 else -0.2
  
  if (type == "var") {
    stopifnot(is.character(select), length(select) == 1)
    v <- select[1]
    stopifnot(v %in% names(data_raw), !is.null(s2p_traces), !is.null(s2skip_traces))
    
    d_e <- s2p_traces    |> dplyr::filter(.data$variable == v)
    d_s <- s2skip_traces |> dplyr::filter(.data$variable == v)
    
    # naive variance from raw data (complete cases)
    x <- as.numeric(data_raw[[v]])
    x <- x[is.finite(x)]
    naive_val <- if (length(x) >= 2) stats::var(x, na.rm = TRUE) else NA_real_
    
    # optional masking (old NA-based masking if you still want it)
    N_est <- round(stats::median(d_e$np * d_e$T, na.rm = TRUE))
    if (mask_pooled) {
      d_e <- d_e |>
        mask_resonances("s2p", N = N_est) |>
        mask_manual_variance_dips("s2p",
                                  hum_n  = var_manual$hum_n,  hum_tol  = var_manual$hum_tol,
                                  soil_n = var_manual$soil_n, soil_tol = var_manual$soil_tol,
                                  co2_n  = var_manual$co2_n,  co2_tol  = var_manual$co2_tol
        )
      y_e <- "s2p_mask"
    } else {
      y_e <- "s2p"
    }
    if (mask_skip) {
      d_s <- d_s |>
        mask_resonances("s2skip", N = N_est) |>
        mask_manual_variance_dips("s2skip",
                                  hum_n  = var_manual$hum_n,  hum_tol  = var_manual$hum_tol,
                                  soil_n = var_manual$soil_n, soil_tol = var_manual$soil_tol,
                                  co2_n  = var_manual$co2_n,  co2_tol  = var_manual$co2_tol
        )
      y_s <- "s2skip_mask"
    } else {
      y_s <- "s2skip"
    }
    
    d_e <- d_e |> dplyr::filter(np > x_limits[1], np < x_limits[2])
    d_s <- d_s |> dplyr::filter(np > x_limits[1], np < x_limits[2])
    
    rects <- build_T_rects_np(
      df = d_e,
      T_col = "T", np_col = "np",
      highlight_Ts          = highlight_Ts,
      highlight_T_halfwidth = highlight_T_halfwidth,
      highlight_labels      = highlight_labels,
      x_limits = x_limits
    )
    
    title_auto <- bquote("Variance — " * .(as.character(.pretty_var(v))))
    
    p <- ggplot() +
      { if (!is.null(rects))
        geom_rect(
          data = rects,
          inherit.aes = FALSE,
          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
          fill  = col_ess,
          alpha = shade_alpha
        )
      } +
      geom_line(data = d_s, aes(x = np, y = .data[[y_s]], color = "Skip", linetype = "Skip"),
                linewidth = lwd_skip, na.rm = TRUE) +
      geom_line(data = d_e, aes(x = np, y = .data[[y_e]], color = "ESS",  linetype = "ESS"),
                linewidth = lwd_ess,  na.rm = TRUE) +
      { if (is.finite(naive_val)) geom_hline(aes(yintercept = naive_val, color = "Naive", linetype = "Naive"),
                                             linewidth = lwd_naive, show.legend = TRUE) } +
      { if (!is.null(rects) && label_T_boxes)
        geom_text(
          data = rects,
          inherit.aes = FALSE,
          aes(x = x_label, label = label),
          y = y_pos,
          vjust = vjust_lab,
          size = 2.5
        )
      } +
    scale_x_sqrt(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_breaks,
      name   = expression(n[t])
    ) +
    scale_color_manual(
      name   = NULL,
      breaks = c("ESS","Skip","Naive"),
      values = c("ESS" = col_ess, "Skip" = col_skip, "Naive" = col_naive)
    ) +
    scale_linetype_manual(
      name   = NULL,
      breaks = c("ESS","Skip","Naive"),
      values = c("ESS" = "solid", "Skip" = "solid", "Naive" = "dashed")
    ) +
    labs(
      x = expression(n[t]),
      y = NULL,
      title = if (is.null(title_prefix)) title_auto else
        bquote(.(title_prefix) ~ " — " ~ .(title_auto))
    ) +
    theme_bw() +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title         = element_text(hjust = 0.5, face = "bold"),
      legend.position    = "top"
    )
    
  } else {  # covariance
    stopifnot(is.character(select), length(select) == 2)
    a <- select[1]; b <- select[2]
    stopifnot(all(c(a,b) %in% names(data_raw)), !is.null(sxy_traces), !is.null(sXYskip_traces))
    
    d_e <- sxy_traces     |> dplyr::filter(.data$pair_key %in% c(paste(a,b,sep="_"), paste(b,a,sep="_")))
    d_s <- sXYskip_traces |> dplyr::filter(.data$pair_key %in% c(paste(a,b,sep="_"), paste(b,a,sep="_")))
    
    # naive covariance from raw data (complete pairs)
    x <- as.numeric(data_raw[[a]]); y <- as.numeric(data_raw[[b]])
    ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
    naive_val <- if (length(x) >= 2) stats::cov(x, y) else NA_real_
    
    # optional masking (old NA-based)
    N_est <- round(stats::median(d_e$np * d_e$T, na.rm = TRUE))
    if (mask_pooled) {
      d_e <- d_e |>
        mask_resonances("sXYp", N = N_est) |>
        mask_manual_cov_dips("sXYp", windows = cov_windows)
      y_e <- "sXYp_mask"
    } else {
      y_e <- "sXYp"
    }
    if (mask_skip) {
      d_s <- d_s |>
        mask_resonances("sXYskip", N = N_est) |>
        mask_manual_cov_dips("sXYskip", windows = cov_windows)
      y_s <- "sXYskip_mask"
    } else {
      y_s <- "sXYskip"
    }
    
    d_e <- d_e |> dplyr::filter(np > x_limits[1], np < x_limits[2])
    d_s <- d_s |> dplyr::filter(np > x_limits[1], np < x_limits[2])
    
    rects <- build_T_rects_np(
      df = d_e,
      T_col = "T", np_col = "np",
      highlight_Ts          = highlight_Ts,
      highlight_T_halfwidth = highlight_T_halfwidth,
      highlight_labels      = highlight_labels,
      x_limits = x_limits
    )
    
    p <- ggplot() +
      { if (!is.null(rects))
        geom_rect(
          data = rects,
          inherit.aes = FALSE,
          aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
          fill  = col_ess,
          alpha = shade_alpha
        )
      } +
      geom_line(data = d_s, aes(x = np, y = .data[[y_s]], color = "Skip", linetype = "Skip"),
                linewidth = lwd_skip, na.rm = TRUE) +
      geom_line(data = d_e, aes(x = np, y = .data[[y_e]], color = "ESS",  linetype = "ESS"),
                linewidth = lwd_ess,  na.rm = TRUE) +
      { if (is.finite(naive_val)) geom_hline(aes(yintercept = naive_val, color = "Naive", linetype = "Naive"),
                                             linewidth = lwd_naive, show.legend = TRUE) } +
      { if (!is.null(rects) && label_T_boxes)
        geom_text(
          data = rects,
          inherit.aes = FALSE,
          aes(x = x_label, label = label),
          y = y_pos,
          vjust = vjust_lab,
          size = 2.5
        )
      } +
      scale_x_sqrt(
        limits = x_limits,
        breaks = x_breaks,
        labels = x_breaks,
        name   = expression(n[t])
      ) +
      scale_color_manual(
        name   = NULL,
        breaks = c("ESS","Skip","Naive"),
        values = c("ESS" = col_ess, "Skip" = col_skip, "Naive" = col_naive)
      ) +
      scale_linetype_manual(
        name   = NULL,
        breaks = c("ESS","Skip","Naive"),
        values = c("ESS" = "solid", "Skip" = "solid", "Naive" = "dashed")
      ) +
      labs(
        x = expression(n[t]),
        y = NULL
      ) +
      theme_bw() +
      theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title         = element_text(hjust = 0.5, face = "bold"),
        legend.position    = "none"
      )
    
  }
  
  if (!is.null(filename)) ggsave(filename, p, width = width, height = height, dpi = dpi)
  p
}

# build skip traces
#s2skip_traces <- compute_s2skip_traces(phase1, vars)
#sXYskip_traces <- compute_pairwise_cov_skip(phase1, vars)

p_humid_soil_naive <- plot_single_ess_skip_naive(
  #select = c("humidity_pct","co2_ppm"),
  select = c("humidity_pct","soil_temp_F"),
  #select = c("co2_ppm","soil_temp_F"),
  data_raw = phase1,
  s2p_traces = s2p_traces, sxy_traces = sxy_traces,
  s2skip_traces = s2skip_traces, sXYskip_traces = sXYskip_traces,
  type = "cov",
  show_T_boxes   = TRUE,
  label_position = "top",   
  filename = "humid_soil_ess_skip_naive.png"
)

#Table of variability in skip
library(dplyr)
library(purrr)
library(tibble)

# ---- helpers reused ----
.ess_var_at_T <- function(x, T) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  N <- length(x); if (N < 2) return(NA_real_)
  vals <- vapply(1:T, function(i){
    idx <- seq.int(i, N, by = T)
    if (length(idx) >= 2) stats::var(x[idx]) else NA_real_
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

.ess_cov_at_T <- function(x, y, T) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n <- length(x); if (n < 2) return(NA_real_)
  vals <- vapply(1:T, function(i){
    idx <- seq.int(i, n, by = T)
    if (length(idx) >= 2) stats::cov(x[idx], y[idx]) else NA_real_
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

.skip_var_all_offsets <- function(x, T) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  N <- length(x)
  vapply(1:T, function(i){
    idx <- seq.int(i, N, by = T)
    if (length(idx) >= 2) stats::var(x[idx]) else NA_real_
  }, numeric(1))
}

.skip_cov_all_offsets <- function(x, y, T) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n <- length(x)
  vapply(1:T, function(i){
    idx <- seq.int(i, n, by = T)
    if (length(idx) >= 2) stats::cov(x[idx], y[idx]) else NA_real_
  }, numeric(1))
}

# ---- main: table over a grid of T values ----
# select: "varname"  OR  c("var_x","var_y")
# data: data.frame containing those columns
# T_seq: optional explicit vector of T's; otherwise seq(T_min, T_max, by = T_by)
# skip_offset: which offset to report for `skip_est` (1..T)
table_over_T <- function(select, data,
                         T_seq = NULL, T_by = 100, T_min = 10, T_max = NULL,
                         skip_offset = 1) {
  stopifnot(is.data.frame(data))
  
  if (length(select) == 1) {
    v <- select[1]; stopifnot(v %in% names(data))
    x <- as.numeric(data[[v]]); x <- x[is.finite(x)]
    N <- length(x); if (is.null(T_max)) T_max <- max(1L, floor(N/10))
    if (is.null(T_seq)) T_seq <- seq(T_min, T_max, by = T_by)
    T_seq <- T_seq[T_seq >= 1]
    naive <- if (length(x) >= 2) stats::var(x) else NA_real_
    
    rows <- purrr::map_dfr(T_seq, function(Tcur){
      if (!is.finite(Tcur) || Tcur < 1) return(tibble::tibble())
      s_all <- .skip_var_all_offsets(x, Tcur)
      off   <- max(1L, min(skip_offset, Tcur))
      idx   <- seq.int(off, N, by = Tcur)
      tibble::tibble(
        type      = "variance",
        variable  = v,
        T         = as.integer(Tcur),
        np        = N / Tcur,
        naive     = naive,
        pooled    = .ess_var_at_T(x, Tcur),
        skip_est  = if (length(idx) >= 2) stats::var(x[idx]) else NA_real_,
        skip_min  = suppressWarnings(min(s_all, na.rm = TRUE)),
        skip_max  = suppressWarnings(max(s_all, na.rm = TRUE)),
        n_used    = length(idx)
      )
    })
    return(rows)
    
  } else if (length(select) == 2) {
    a <- select[1]; b <- select[2]
    stopifnot(all(c(a,b) %in% names(data)))
    x <- as.numeric(data[[a]]); y <- as.numeric(data[[b]])
    ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
    n <- length(x); if (is.null(T_max)) T_max <- max(1L, floor(n/10))
    if (is.null(T_seq)) T_seq <- seq(T_min, T_max, by = T_by)
    T_seq <- T_seq[T_seq >= 1]
    naive <- if (length(x) >= 2) stats::cov(x, y) else NA_real_
    
    rows <- purrr::map_dfr(T_seq, function(Tcur){
      if (!is.finite(Tcur) || Tcur < 1) return(tibble::tibble())
      s_all <- .skip_cov_all_offsets(x, y, Tcur)
      off   <- max(1L, min(skip_offset, Tcur))
      idx   <- seq.int(off, n, by = Tcur)
      tibble::tibble(
        type      = "covariance",
        var_x     = a,
        var_y     = b,
        T         = as.integer(Tcur),
        np        = n / Tcur,
        naive     = naive,
        pooled    = .ess_cov_at_T(x, y, Tcur),
        skip_est  = if (length(idx) >= 2) stats::cov(x[idx], y[idx]) else NA_real_,
        skip_min  = suppressWarnings(min(s_all, na.rm = TRUE)),
        skip_max  = suppressWarnings(max(s_all, na.rm = TRUE)),
        n_used    = length(idx)
      )
    })
    return(rows)
    
  } else {
    stop("`select` must be a single variable (variance) or length-2 vector (covariance).")
  }
}

library(xtable)
library(dplyr)
print_table_over_T_xtable <- function(tab,
                                      caption = NULL,
                                      label   = "tab:overT",
                                      digits_T   = 0,
                                      digits_n   = 0,
                                      digits_est = 4,
                                      use_booktabs = TRUE) {
  stopifnot(is.data.frame(tab))
  
  # figure out naive estimate (assume constant within tab)
  naive_vals <- unique(tab$naive[is.finite(tab$naive)])
  naive_val  <- if (length(naive_vals)) naive_vals[1] else NA_real_
  
  # type for caption: variance / covariance
  type_vals <- unique(tab$type)
  kind <- if (length(type_vals) == 1L) {
    if (type_vals == "variance") "Variance"
    else if (type_vals == "covariance") "Covariance"
    else "Summary"
  } else {
    "Summary"
  }
  
  if (is.null(caption)) {
    if (is.finite(naive_val)) {
      caption <- sprintf(
        "%s across $T$ with ESS and skip-based estimates. Naive estimate: %.4f.",
        kind, naive_val
      )
    } else {
      caption <- sprintf(
        "%s across $T$ with ESS and skip-based estimates.", kind
      )
    }
  }
  
  # Build the body: exactly the columns you want
  df <- tab %>%
    dplyr::transmute(
      T        = as.integer(T),
      n        = as.integer(n_used),  # labeled as n
      ESS      = pooled,
      Skip_min = skip_min,
      Skip_est = skip_est,
      Skip_max = skip_max
    )
  
  # xtable object
  xt <- xtable::xtable(
    df,
    caption = caption,
    label   = label,
    align   = c("r","r","r","r","r","r","r")  # 6 numeric cols
  )
  
  # digits: rownames, T, n, ESS, min, est, max
  digits_vec <- c(
    0,            # rownames
    digits_T,     # T
    digits_n,     # n
    digits_est,   # ESS
    digits_est,   # Skip_min
    digits_est,   # Skip_est
    digits_est    # Skip_max
  )
  attr(xt, "digits") <- digits_vec
  
  # multicolumn header:
  #   T & n & ESS & \multicolumn{3}{c}{Skip-Sampling} \\
  #   \cline{4-6}
  #    &   &     & Min & Offset=1 Estimate & Max \\
  if (use_booktabs) {
    addtorow <- list(
      pos = list(0, nrow(df)),
      command = c(
        paste0(
          "\\toprule\n",
          "T & n & ESS & \\multicolumn{3}{c}{Skip-Sampling} \\\\\n",
          "\\cline{4-6}\n",
          " &  &  & Min & Offset=1 Estimate & Max \\\\\n",
          "\\midrule\n"
        ),
        "\\bottomrule\n"
      )
    )
  } else {
    addtorow <- list(
      pos = list(0),
      command = c(
        paste0(
          "T & n & ESS & \\multicolumn{3}{c}{Skip-Sampling} \\\\\n",
          "\\cline{4-6}\n",
          " &  &  & Min & Offset=1 Estimate & Max \\\\\n",
          "\\hline\n"
        )
      )
    )
  }
  
  print(
    xt,
    include.rownames = FALSE,
    include.colnames = FALSE,          # we’re supplying our own header
    sanitize.text.function = identity, # keep LaTeX as-is
    floating = TRUE,
    caption.placement = "top",
    hline.after = NULL,                # handled in add.to.row
    add.to.row = addtorow,
    booktabs = use_booktabs
  )
}


tab_cov_hum_soil <- table_over_T(
  select = c("humidity_pct", "soil_temp_F"),
  data   = phase1,
  T_seq  = c(10, 50, 100, 200, 500, 1000),
  skip_offset = 1
)

print_table_over_T_xtable(
  tab_cov_hum_soil,
  label   = "tab:cov_hum_soil",
  digits_T = 0,
  digits_n = 0,
  digits_est = 2
)

tab_cov_hum_co2 <- table_over_T(
  select = c("humidity_pct", "co2_ppm"),
  data   = phase1,
  T_seq  = c(10, 50, 100, 200, 500, 1000),
  skip_offset = 1
)

tab_cov_co2_soil <- table_over_T(
  select = c("co2_ppm", "soil_temp_F"),
  data   = phase1,
  T_seq  = c(10, 50, 100, 200, 500, 1000),
  skip_offset = 1
)


tab_var <- table_over_T(
  select = "co2_ppm",
  data   = phase1,
  T_seq  = c(10, 50, 100, 200, 300),
  skip_offset = 1
)
print(tab_var)

# install.packages("xtable") # if needed
library(dplyr)
library(xtable)

# ---- format & print with xtable ----
print_table_over_T_xtable <- function(tab,
                                      caption = NULL,
                                      label = "tab:overT",
                                      digits_np = 1,   # digits for n = N/T
                                      digits_est = 4,  # digits for estimates
                                      use_booktabs = TRUE) {
  stopifnot(is.data.frame(tab))
  
  # If it's covariance, you likely have var_x/var_y columns; we don’t need them in the table body.
  has_cov <- all(c("var_x","var_y") %in% names(tab))
  var_title <- if (has_cov) paste0(tab$var_x[1], " vs ", tab$var_y[1]) else tab$variable[1]
  
  # Build the table body
  df <- tab %>%
    transmute(
      T                    = as.integer(T),
      n                    = np,          # N/T
      naive                = naive,
      pooled               = pooled,      # ESS (pooled)
      skip                 = skip_est,    # skip for the chosen offset
      skip_min             = pmin(skip_min, skip_max, na.rm = TRUE),
      skip_max             = pmax(skip_min, skip_max, na.rm = TRUE),
      n_used_ref           = as.integer(n_used_ref)
    )
  
  # Column names (math-ready)
  colnames(df) <- c(
    "$T$",
    "$n$",
    "\\textit{Naive}",
    "\\textit{ESS (pooled)}",
    "\\textit{Skip (offset 1)}",
    "\\min(\\textit{Skip})",
    "\\max(\\textit{Skip})",
    "$n_{\\text{used}}$"
  )
  
  # xtable object
  if (is.null(caption)) {
    cap_kind <- if (has_cov) "Covariance" else "Variance"
    caption <- paste0(cap_kind, " across $T$ for ", var_title,
                      " (Naive, ESS pooled, Skip (offset 1), and Skip range).")
  }
  
  xt <- xtable(df,
               caption = caption,
               label   = label,
               align   = c("l", "r","r","r","r","r","r","r","r"))
  
  # Digits: first entry is for rownames (0), then per column
  digits_vec <- c(
    0,                 # rownames
    0,                 # T
    digits_np,         # n
    digits_est,        # naive
    digits_est,        # pooled
    digits_est,        # skip
    digits_est,        # skip_min
    digits_est,        # skip_max
    0                  # n_used_ref
  )
  attr(xt, "digits") <- digits_vec
  
  # Printing
  if (use_booktabs) {
    addtorow <- list(
      pos = list(0, nrow(df)),
      command = c("\\toprule\n", "\\bottomrule\n")
    )
    print(xt,
          include.rownames = FALSE,
          include.colnames = TRUE,
          sanitize.text.function = identity,      # keep math in headers
          floating = TRUE, caption.placement = "top",
          hline.after = c(),                      # let booktabs handle rules
          add.to.row = addtorow)
  } else {
    print(xt,
          include.rownames = FALSE,
          include.colnames = TRUE,
          sanitize.text.function = identity,
          floating = TRUE, caption.placement = "top",
          hline.after = c(-1, 0, nrow(df)))       # standard hlines
  }
}

# ---- example using your table ----
print_table_over_T_xtable(
  tab_cov_hum_soil,
  caption = "Covariance across $T$ for Humidity versus Soil Temperature.",
  label   = "tab:cov_hum_soil",
  digits_np = 1,
  digits_est = 1,
  use_booktabs = TRUE
)

print_table_over_T_xtable(
  tab_cov_hum_co2,
  caption = "Covariance across $T$ for Humidity versus CO$_2$.",
  label   = "tab:cov_hum_co2",
  digits_np = 1,
  digits_est = 1,
  use_booktabs = TRUE
)

print_table_over_T_xtable(
  tab_cov_co2_soil,
  caption = "Covariance across $T$ for CO$_2$ versus Soil Temperature.",
  label   = "tab:cov_co2_soil",
  digits_np = 1,
  digits_est = 1,
  use_booktabs = TRUE
)



# Find Ts and estimates

# --- deps
library(dplyr)
library(tidyr)
library(purrr)
library(zoo)

# ---------- ACF/CCF gate identical to your sim code ----------
first_below_thresh_lag <- function(x, y = NULL, lag_max = 40) {
  ac <- if (is.null(y)) {
    stats::acf(x, lag.max = lag_max, plot = FALSE, na.action = stats::na.pass)
  } else {
    stats::ccf(x, y, lag.max = lag_max, plot = FALSE, na.action = stats::na.pass)
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
  rr <- range(x, na.rm = TRUE)
  diff(rr)
}

# pick the nearest non-NA value to a target T
.value_at_T_or_nearest <- function(df, ycol, T_target) {
  y <- df[[ycol]]
  if (!is.na(match(T_target, df$T)) && !is.na(y[df$T == T_target][1])) {
    i <- which(df$T == T_target)[1]
    return(list(T = df$T[i], np = df$np[i], value = y[i]))
  }
  ok <- which(!is.na(y))
  if (!length(ok)) return(list(T = NA_integer_, np = NA_real_, value = NA_real_))
  i <- ok[ which.min(abs(df$T[ok] - T_target)) ]
  list(T = df$T[i], np = df$np[i], value = y[i])
}

# choose T from one pooled series (already masked & filtered)
.chooseT_one <- function(df, ycol_mask,
                         method = c("acf","acf_roll","roll"),
                         window = 20, frac = 0.30,
                         lag_max = 40,
                         X = NULL, Y = NULL) {
  method <- match.arg(method)
  
  # masked series to scan
  series <- df[[ycol_mask]]
  n      <- length(series)
  
  if (n == 0L || all(is.na(series))) {
    return(list(T_roll = NA_integer_, T_corr = NA_integer_, T_star = NA_integer_))
  }
  
  w <- min(window, n)
  
  # total range and threshold
  total_range <- diff(range(series, na.rm = TRUE))
  thr         <- frac * total_range
  
  # rolling range over indices 1..n (window width w)
  roll_rng <- zoo::rollapplyr(
    series, width = w,
    FUN = .roll_range, partial = FALSE, fill = NA_real_
  )
  
  exceed_idx <- which(roll_rng > thr)
  
  # ---------- pick index for T_roll using the "previous window median" rule ----------
  if (length(exceed_idx)) {
    # first window whose range exceeds threshold
    i_pick  <- exceed_idx[1]               # right edge of triggering window
    j_start <- i_pick - w + 1              # left edge of triggering window
    j_end   <- i_pick
    
    # window that ABUTS the triggering one: [prev_start, prev_end]
    prev_end <- j_start - 1
    if (prev_end >= 1L) {
      prev_start <- max(1L, prev_end - w + 1L)
      idx_win    <- prev_start:prev_end
    } else {
      # no full previous window; default to the earliest data we have
      idx_win <- 1:j_end
    }
  } else {
    # no alarm: use the last window and take its median point
    if (n <= w) {
      idx_win <- 1:n
    } else {
      idx_win <- (n - w + 1L):n
    }
  }
  
  window_y <- series[idx_win]
  
  if (all(is.na(window_y))) {
    # degenerate case: fall back to the first index of the window
    i_roll <- idx_win[1]
  } else {
    med <- stats::median(window_y, na.rm = TRUE)
    rel <- which(!is.na(window_y))
    k   <- rel[ which.min(abs(window_y[rel] - med)) ]  # position within the window
    i_roll <- idx_win[k]                               # global index in df
  }
  
  T_roll <- df$T[i_roll]
  
  # ---------- ACF/CCF gate (unchanged) ----------
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


# Drop any rows whose T lies in one of the T0 ± halfwidth windows
drop_T_windows <- function(df,
                           T_col = "T",
                           highlight_Ts,
                           highlight_T_halfwidth) {
  if (is.null(highlight_Ts) || !length(highlight_Ts)) return(df)
  
  # Recycle halfwidth if scalar
  if (length(highlight_T_halfwidth) == 1L &&
      length(highlight_Ts) > 1L) {
    highlight_T_halfwidth <- rep(highlight_T_halfwidth, length(highlight_Ts))
  }
  stopifnot(length(highlight_Ts) == length(highlight_T_halfwidth))
  
  Tvals <- df[[T_col]]
  bad <- rep(FALSE, length(Tvals))
  
  for (k in seq_along(highlight_Ts)) {
    T0 <- highlight_Ts[k]
    hw <- highlight_T_halfwidth[k]
    if (!is.finite(T0) || !is.finite(hw)) next
    bad <- bad | (Tvals >= (T0 - hw) & Tvals <= (T0 + hw))
  }
  
  df[!bad, , drop = FALSE]
}


# =================== MAIN: choose Ts from greenhouse traces ===================
choose_greenhouse_Ts <- function(s2p_traces, sxy_traces, data_raw = NULL,
                                 methods = "roll",
                                 chooseT_window = 20, chooseT_frac = 0.30,
                                 lag_max = NULL,
                                 x_limits = c(10, 300),
                                 # NEW: T windows to exclude entirely
                                 highlight_Ts          = c(480, 720, 960, 1440, 2160, 2880, 3600, 4320),
                                 highlight_T_halfwidth = c(20,  25,  20,  95,   40,   175,  60,   60)) {
  
  stopifnot(all(c("T","np","s2p","variable") %in% names(s2p_traces)))
  stopifnot(all(c("T","np","sXYp","var_x","var_y","pair_key") %in% names(sxy_traces)))
  methods <- intersect(unique(methods), c("acf","acf_roll","roll"))
  if (!length(methods)) stop("methods must include at least one of: 'acf','acf_roll','roll'.")
  
  if (is.null(lag_max)) lag_max <- 40L
  
  # ----- VARIANCES: filter in n-space, then drop T-windows, then define s2p_mask -----
  vdat <- s2p_traces %>%
    dplyr::filter(np > x_limits[1], np < x_limits[2]) %>%
    drop_T_windows(T_col = "T",
                   highlight_Ts = highlight_Ts,
                   highlight_T_halfwidth = highlight_T_halfwidth) %>%
    dplyr::mutate(s2p_mask = s2p)   # mask column used by .chooseT_one
  
  # ----- COVARIANCES: same idea -----
  cdat <- sxy_traces %>%
    dplyr::filter(np > x_limits[1], np < x_limits[2]) %>%
    drop_T_windows(T_col = "T",
                   highlight_Ts = highlight_Ts,
                   highlight_T_halfwidth = highlight_T_halfwidth) %>%
    dplyr::mutate(sXYp_mask = sXYp)
  
  # ---------- choose Ts per variable ----------
  var_choices <- purrr::map_dfr(unique(vdat$variable), function(v) {
    dfv <- vdat %>%
      dplyr::filter(variable == v) %>%
      dplyr::arrange(T)
    
    if (!nrow(dfv)) {
      return(tibble::tibble(
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
    
    purrr::map_dfr(methods, function(m) {
      choice <- .chooseT_one(dfv, "s2p_mask", method = m,
                             window = chooseT_window, frac = chooseT_frac,
                             lag_max = lag_max, X = X, Y = NULL)
      
      picked <- .value_at_T_or_nearest(dfv, "s2p_mask", choice$T_star)
      
      naive_row <- s2p_traces %>% dplyr::filter(variable == v, T == 1)
      naive_val <- if (nrow(naive_row)) naive_row$s2p[1] else NA_real_
      
      tibble::tibble(
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
  
  # ---------- choose Ts per pair ----------
  cov_choices <- purrr::map_dfr(unique(cdat$pair_key), function(pk) {
    dfp <- cdat %>%
      dplyr::filter(pair_key == pk) %>%
      dplyr::arrange(T)
    
    if (!nrow(dfp)) {
      # recover var names from original sxy_traces for a consistent row
      tmp <- sxy_traces %>% dplyr::filter(pair_key == pk) %>% dplyr::slice(1)
      vx <- tmp$var_x[1]; vy <- tmp$var_y[1]
      return(tibble::tibble(
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
    
    vx <- dfp$var_x[1]; vy <- dfp$var_y[1]
    
    # raw series for CCF gate
    X <- if (!is.null(data_raw)) as.numeric(data_raw[[vx]]) else NULL
    Y <- if (!is.null(data_raw)) as.numeric(data_raw[[vy]]) else NULL
    if (!is.null(X) && !is.null(Y)) {
      ok <- is.finite(X) & is.finite(Y)
      X <- X[ok]; Y <- Y[ok]
    } else {
      X <- Y <- NULL
    }
    
    purrr::map_dfr(methods, function(m) {
      choice <- .chooseT_one(dfp, "sXYp_mask", method = m,
                             window = chooseT_window, frac = chooseT_frac,
                             lag_max = lag_max, X = X, Y = Y)
      
      picked <- .value_at_T_or_nearest(dfp, "sXYp_mask", choice$T_star)
      
      naive_row <- sxy_traces %>% dplyr::filter(pair_key == pk, T == 1)
      naive_val <- if (nrow(naive_row)) naive_row$sXYp[1] else NA_real_
      
      tibble::tibble(
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
  
  list(var_choices = var_choices, cov_choices = cov_choices)
}

# Optional: pair-specific manual masks for covariance (edit as you like)
highlight_Ts          <- c(480, 720, 960, 1440, 2160, 2880, 3600, 4320)
highlight_T_halfwidth <- c(20,  25,  20,   95,   40,   175,   60,   60)

res_T <- choose_greenhouse_Ts(
  s2p_traces, sxy_traces, data_raw = phase1,
  methods = c("roll"),
  chooseT_window = 30, chooseT_frac = 0.2,  x_limits = c(10, 300),
  highlight_Ts          = highlight_Ts,
  highlight_T_halfwidth = highlight_T_halfwidth
)

res_T$var_choices
res_T$cov_choices

plot_pooled_var_cov(
  s2p_traces,
  sxy_traces,
  highlight_Ts          = highlight_Ts,
  highlight_T_halfwidth = highlight_T_halfwidth,
  res_T      = res_T,
  sel_method = "roll",         
  filename   = "greenhouse_var_cov_with_Tdots.png"
)


