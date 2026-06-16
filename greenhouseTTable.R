ess_var_at_T <- function(x, T) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  N <- length(x)
  if (N < 2) return(NA_real_)
  
  vals <- vapply(seq_len(T), function(i) {
    idx <- seq.int(i, N, by = T)
    if (length(idx) >= 2) stats::var(x[idx]) else NA_real_
  }, numeric(1))
  
  mean(vals, na.rm = TRUE)
}

ess_cov_at_T <- function(x, y, T) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x  <- x[ok];       y <- y[ok]
  n  <- length(x)
  if (n < 2) return(NA_real_)
  
  vals <- vapply(seq_len(T), function(i) {
    idx <- seq.int(i, n, by = T)
    if (length(idx) >= 2) stats::cov(x[idx], y[idx]) else NA_real_
  }, numeric(1))
  
  mean(vals, na.rm = TRUE)
}

skip_var_all_offsets <- function(x, T) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  N <- length(x)
  
  vapply(seq_len(T), function(i) {
    idx <- seq.int(i, N, by = T)
    if (length(idx) >= 2) stats::var(x[idx]) else NA_real_
  }, numeric(1))
}

skip_cov_all_offsets <- function(x, y, T) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x  <- x[ok];       y <- y[ok]
  n  <- length(x)
  
  vapply(seq_len(T), function(i) {
    idx <- seq.int(i, n, by = T)
    if (length(idx) >= 2) stats::cov(x[idx], y[idx]) else NA_real_
  }, numeric(1))
}

# table over a grid of T values 
# select:
#   - "varname"             -> variance table
#   - c("var_x", "var_y")   -> covariance table
#
# data: data.frame containing those columns
# T_seq: explicit vector of T's, or built from (T_min, T_max, T_by)
# skip_offset: which offset (1..T) to report as "Skip" in the table
table_over_T <- function(select,
                         data,
                         T_seq      = NULL,
                         T_by       = 100,
                         T_min      = 10,
                         T_max      = NULL,
                         skip_offset = 1) {
  stopifnot(is.data.frame(data))
  
  # variance case 
  if (length(select) == 1) {
    v <- select[1]
    stopifnot(v %in% names(data))
    
    x <- as.numeric(data[[v]])
    x <- x[is.finite(x)]
    N <- length(x)
    if (N < 2) stop("Not enough data for variance.")
    
    if (is.null(T_max)) T_max <- max(1L, floor(N / 10))
    if (is.null(T_seq)) T_seq <- seq(T_min, T_max, by = T_by)
    T_seq <- as.integer(T_seq[T_seq >= 1 & T_seq <= T_max])
    if (!length(T_seq)) stop("No valid T values in T_seq.")
    
    naive <- stats::var(x)
    
    rows <- purrr::map_dfr(T_seq, function(Tcur) {
      s_all <- skip_var_all_offsets(x, Tcur)
      off   <- max(1L, min(skip_offset, Tcur))
      idx   <- seq.int(off, N, by = Tcur)
      
      tibble::tibble(
        type      = "variance",
        variable  = v,
        T         = as.integer(Tcur),
        np        = N / Tcur,
        naive     = naive,
        pooled    = ess_var_at_T(x, Tcur),
        skip_est  = if (length(idx) >= 2) stats::var(x[idx]) else NA_real_,
        skip_min  = suppressWarnings(min(s_all, na.rm = TRUE)),
        skip_max  = suppressWarnings(max(s_all, na.rm = TRUE)),
        n_used    = length(idx)
      )
    })
    
    return(rows)
  }
  
  # covariance case
  if (length(select) == 2) {
    a <- select[1]
    b <- select[2]
    stopifnot(all(c(a, b) %in% names(data)))
    
    x <- as.numeric(data[[a]])
    y <- as.numeric(data[[b]])
    ok <- is.finite(x) & is.finite(y)
    x  <- x[ok];       y <- y[ok]
    n  <- length(x)
    if (n < 2) stop("Not enough data for covariance.")
    
    if (is.null(T_max)) T_max <- max(1L, floor(n / 10))
    if (is.null(T_seq)) T_seq <- seq(T_min, T_max, by = T_by)
    T_seq <- as.integer(T_seq[T_seq >= 1 & T_seq <= T_max])
    if (!length(T_seq)) stop("No valid T values in T_seq.")
    
    naive <- stats::cov(x, y)
    
    rows <- purrr::map_dfr(T_seq, function(Tcur) {
      s_all <- skip_cov_all_offsets(x, y, Tcur)
      off   <- max(1L, min(skip_offset, Tcur))
      idx   <- seq.int(off, n, by = Tcur)
      
      tibble::tibble(
        type      = "covariance",
        var_x     = a,
        var_y     = b,
        T         = as.integer(Tcur),
        np        = n / Tcur,
        naive     = naive,
        pooled    = ess_cov_at_T(x, y, Tcur),
        skip_est  = if (length(idx) >= 2) stats::cov(x[idx], y[idx]) else NA_real_,
        skip_min  = suppressWarnings(min(s_all, na.rm = TRUE)),
        skip_max  = suppressWarnings(max(s_all, na.rm = TRUE)),
        n_used    = length(idx)
      )
    })
    
    return(rows)
  }
  
  stop("`select` must be a single name (variance) or length-2 (covariance).")
}

print_table_over_T_xtable <- function(tab,
                                      caption    = NULL,
                                      label      = "tab:overT",
                                      digits_T   = 0,
                                      digits_n   = 0,
                                      digits_est = 2,
                                      use_booktabs = TRUE) {
  stopifnot(is.data.frame(tab))
  
  # Figure out if this is var or cov for default caption
  type_vals <- unique(tab$type)
  has_cov   <- all(c("var_x", "var_y") %in% names(tab))
  
  var_title <- if (has_cov) {
    paste0(tab$var_x[1], " vs ", tab$var_y[1])
  } else {
    tab$variable[1]
  }
  
  if (is.null(caption)) {
    cap_kind <- if (length(type_vals) == 1L && type_vals == "covariance") {
      "Covariance"
    } else if (length(type_vals) == 1L && type_vals == "variance") {
      "Variance"
    } else {
      "Summary"
    }
    caption <- paste0(
      cap_kind, " across $T$ for ", var_title, "."
    )
  }
  
  # Build exactly the 6 columns you want in the body:
  # T, n (=N/T), ESS (=pooled), Skip min, Skip (offset 1), Skip max
  df <- tab %>%
    dplyr::transmute(
      T        = as.integer(T),
      n        = round(np),        # N / T, rounded like your example (e.g. 4158)
      ESS      = pooled,
      Skip_min = skip_min,
      Skip_est = skip_est,
      Skip_max = skip_max
    )
  
  # Create xtable object with 6 numeric columns
  xt <- xtable::xtable(
    df,
    caption = caption,
    label   = label,
    align   = c("r", "r", "r", "r", "r", "r", "r")  # 6 cols -> length = 7
  )
  
  # digits: rownames, T, n, ESS, Min, Offset1, Max
  digits_vec <- c(
    0,           # rownames
    digits_T,    # T
    digits_n,    # n
    digits_est,  # ESS
    digits_est,  # Min
    digits_est,  # Offset 1
    digits_est   # Max
  )
  attr(xt, "digits") <- digits_vec
  
  if (use_booktabs) {
    # Header + footer exactly in your style
    addtorow <- list(
      pos = list(0, nrow(df)),
      command = c(
        paste0(
          "\\toprule\n",
          "$T$ & $n_t$ & ESS & \\multicolumn{3}{c}{Skip-Sampling} \\\\\n",
          "\\cline{4-6}\n",
          " &  &  & Min & Offset 1 & Max \\\\\n",
          "\\midrule\n"
        ),
        "\\bottomrule\n"
      )
    )
    
    print(
      xt,
      include.rownames       = FALSE,
      include.colnames       = FALSE,      # we provide our own header
      sanitize.text.function = identity,   # keep $T$, etc.
      floating               = TRUE,
      caption.placement      = "top",
      hline.after            = NULL,
      add.to.row             = addtorow,
      booktabs               = TRUE
    )
  } else {
    # non-booktabs fallback (simpler header)
    addtorow <- list(
      pos = list(0),
      command = c(
        paste0(
          "$T$ & $n_t$ & ESS & \\multicolumn{3}{c}{Skip-Sampling} \\\\\n",
          "\\cline{4-6}\n",
          " &  &  & Min & Offset 1 & Max \\\\\n",
          "\\hline\n"
        )
      )
    )
    
    print(
      xt,
      include.rownames       = FALSE,
      include.colnames       = FALSE,
      sanitize.text.function = identity,
      floating               = TRUE,
      caption.placement      = "top",
      hline.after            = NULL,
      add.to.row             = addtorow
    )
  }
}


# covariance tables
tab_cov_hum_soil <- table_over_T(
  select      = c("humidity_pct", "soil_temp_F"),
  data        = phase1,
  T_seq       = c(1, 10, 50, 100, 200, 500, 1000),
  skip_offset = 1
)

tab_cov_hum_co2 <- table_over_T(
  select      = c("humidity_pct", "co2_ppm"),
  data        = phase1,
  T_seq       = c(1, 10, 50, 100, 500, 1000),
  skip_offset = 1
)

tab_cov_co2_soil <- table_over_T(
  select      = c("co2_ppm", "soil_temp_F"),
  data        = phase1,
  T_seq       = c(1, 10, 50, 100, 500, 1000),
  skip_offset = 1
)

# variance example
tab_var_co2 <- table_over_T(
  select      = "co2_ppm",
  data        = phase1,
  T_seq       = c(1, 10, 50, 100, 200, 300),
  skip_offset = 1
)

# print as LaTeX tables
# this is table 1
print_table_over_T_xtable(
  tab_cov_hum_soil,
  caption    = "Estimates of Covariance between Humidity and Soil Temperature as $T$ increases.",
  label      = "tab:cov_hum_soil",
  digits_T  = 0,
  digits_n = 0,
  digits_est = 2
)

print_table_over_T_xtable(
  tab_cov_hum_co2,
  caption    = "Estimates of Covariance across $T$ for Humidity versus CO$_2$.",
  label      = "tab:cov_hum_co2",
  digits_n  = 1,
  digits_est = 2
)

print_table_over_T_xtable(
  tab_cov_co2_soil,
  caption    = "Estimates of Covariance across $T$ for CO$_2$ versus Soil Temperature.",
  label      = "tab:cov_co2_soil",
  digits_n  = 1,
  digits_est = 2
)

print_table_over_T_xtable(
  tab_var_co2,
  caption    = "Estimates of Variance across $T$ for CO$_2$.",
  label      = "tab:var_co2",
  digits_n  = 1,
  digits_est = 2
)
