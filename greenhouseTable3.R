print_greenhouse_T_table <- function(res_T,
                                     var_order = c("co2_ppm", "humidity_pct", "soil_temp_F"),
                                     method = "roll") {
  
  vc <- res_T$var_choices
  cc <- res_T$cov_choices
  
  vc <- vc[vc$method == method, , drop = FALSE]
  cc <- cc[cc$method == method, , drop = FALSE]
  
  # Labels for table entries
  var_labels <- c(
    co2_ppm       = "\\text{CO}_2",
    humidity_pct  = "\\text{Humidity}",
    soil_temp_F   = "\\text{Soil temp}"
  )
  
  # Variance rows
  var_rows <- data.frame(
    Entry = paste0("$\\mathrm{Var}(", var_labels[var_order], ")$"),
    T_selected = NA_integer_,
    n_t = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(var_order)) {
    v <- var_order[i]
    row_v <- vc[vc$variable == v, , drop = FALSE]
    
    if (!nrow(row_v)) {
      stop("No var_choices row for variable ", v)
    }
    
    var_rows$T_selected[i] <- row_v$T_selected[1]
    var_rows$n_t[i] <- as.integer(round(row_v$np_selected[1]))
  }
  
  # Covariance rows
  cov_pairs <- list(
    c("co2_ppm", "humidity_pct"),
    c("co2_ppm", "soil_temp_F"),
    c("humidity_pct", "soil_temp_F")
  )
  
  cov_rows <- data.frame(
    Entry = character(length(cov_pairs)),
    T_selected = NA_integer_,
    n_t = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  for (k in seq_along(cov_pairs)) {
    vx <- cov_pairs[[k]][1]
    vy <- cov_pairs[[k]][2]
    
    row_c <- cc[
      (cc$var_x == vx & cc$var_y == vy) |
        (cc$var_x == vy & cc$var_y == vx),
      ,
      drop = FALSE
    ]
    
    if (!nrow(row_c)) {
      stop("No cov_choices row for pair ", vx, " / ", vy)
    }
    
    cov_rows$Entry[k] <- paste0(
      "$\\mathrm{Cov}(",
      var_labels[vx], ",",
      var_labels[vy],
      ")$"
    )
    cov_rows$T_selected[k] <- row_c$T_selected[1]
    cov_rows$n_t[k] <- as.integer(round(row_c$np_selected[1]))
  }
  
  tab <- rbind(var_rows, cov_rows)
  
  # Print LaTeX table
  cat("\\begin{table}[h]\n")
  cat("\\centering\n")
  cat("\\caption{Selected spacings $T$ and resulting within-subset sizes $n_t$ used for the ESS estimates in Table~\\ref{tab:Sigma-greenhouse}.}\n")
  cat("\\label{tab:Sigma-greenhouse-Ts}\n")
  cat("\\renewcommand{\\arraystretch}{1.1}\n")
  cat("\\begin{tabular}{lcc}\n")
  cat("\\toprule\n")
  cat("Entry & $T_{\\text{selected}}$ & $n_t$ \\\\\n")
  cat("\\midrule\n")
  
  for (i in seq_len(nrow(tab))) {
    cat(
      tab$Entry[i], " & ",
      tab$T_selected[i], " & ",
      tab$n_t[i], " \\\\\n",
      sep = ""
    )
  }
  
  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
}

print_greenhouse_T_table(
  res_T,
  var_order = c("co2_ppm", "humidity_pct", "soil_temp_F"),
  method = "roll"
)
