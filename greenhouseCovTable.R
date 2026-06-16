# create table 2: covariance estimates

## build ESS & naive covariance matrices from res_T 

build_cov_matrices_from_resT <- function(res_T,
                                         var_order = c("co2_ppm", "humidity_pct", "soil_temp_F"),
                                         method    = "roll") {
  stopifnot(is.list(res_T),
            all(c("var_choices", "cov_choices") %in% names(res_T)))
  
  vc <- res_T$var_choices
  cc <- res_T$cov_choices
  
  if (!all(c("variable", "method", "estimate", "naive_est") %in% names(vc))) {
    stop("res_T$var_choices does not have expected columns.")
  }
  if (!all(c("var_x", "var_y", "method", "estimate", "naive_est") %in% names(cc))) {
    stop("res_T$cov_choices does not have expected columns.")
  }
  
  vc <- vc[vc$method == method, , drop = FALSE]
  cc <- cc[cc$method == method, , drop = FALSE]
  
  p <- length(var_order)
  Sigma_ess   <- matrix(NA_real_, nrow = p, ncol = p,
                        dimnames = list(var_order, var_order))
  Sigma_naive <- matrix(NA_real_, nrow = p, ncol = p,
                        dimnames = list(var_order, var_order))
  
  ## diagonals: variances
  for (i in seq_len(p)) {
    v  <- var_order[i]
    row_v <- vc[vc$variable == v, , drop = FALSE]
    if (!nrow(row_v)) stop("No var_choices row for variable ", v)
    Sigma_ess[i, i]   <- row_v$estimate[1]
    Sigma_naive[i, i] <- row_v$naive_est[1]
  }
  
  ## off-diagonals: covariances (symmetrized)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (j <= i) next
      vx <- var_order[i]
      vy <- var_order[j]
      row_c <- cc[
        (cc$var_x == vx & cc$var_y == vy) |
          (cc$var_x == vy & cc$var_y == vx),
        ,
        drop = FALSE
      ]
      if (!nrow(row_c)) stop("No cov_choices row for pair ", vx, " / ", vy)
      Sigma_ess[i, j]   <- Sigma_ess[j, i]   <- row_c$estimate[1]
      Sigma_naive[i, j] <- Sigma_naive[j, i] <- row_c$naive_est[1]
    }
  }
  
  list(Sigma_ess = Sigma_ess, Sigma_naive = Sigma_naive)
}

## convert a numeric matrix to a LaTeX bmatrix string 

matrix_to_bmatrix <- function(M, digits = 1) {
  fmt <- function(x) formatC(x, format = "f", digits = digits)
  rows <- apply(M, 1, function(r) paste(fmt(r), collapse = " & "))
  body <- paste0(rows, "\\\\", collapse = "\n")
  paste0("\\begin{bmatrix}\n", body, "\n\\end{bmatrix}")
}

############################
# Get Newey West estimator
# load data
library(readr)
library(sandwich)

greenhouse2 <- read_csv("GreenhouseEnvironmentTimeSeries.csv")

var_order <- c("co2_ppm", "humidity_pct", "soil_temp_F")

# Make sure Newey-West uses the same variable order as the table
greenhouseNW <- as.matrix(greenhouse2[, var_order])

Sigma_NW <- NeweyWest(lm(greenhouseNW ~ 1)) * nrow(greenhouseNW)

# If needed, force names/order to match the table
dimnames(Sigma_NW) <- list(var_order, var_order)
############################


## LaTeX table for Sigma (ESS vs naive) 
print_greenhouse_Sigma_table <- function(res_T,
                                         Sigma_NW,
                                         var_order = c("co2_ppm", "humidity_pct", "soil_temp_F"),
                                         method    = "roll",
                                         digits    = 1) {
  mats <- build_cov_matrices_from_resT(res_T, var_order = var_order, method = method)

  Sigma_ess   <- mats$Sigma_ess
  Sigma_naive <- mats$Sigma_naive

  Sigma_NW <- as.matrix(Sigma_NW)

  if (!all(dim(Sigma_NW) == c(length(var_order), length(var_order)))) {
    stop("Sigma_NW does not have dimensions matching var_order.")
  }

  # Reorder if row/column names are present
  if (!is.null(rownames(Sigma_NW)) && !is.null(colnames(Sigma_NW)) &&
      all(var_order %in% rownames(Sigma_NW)) &&
      all(var_order %in% colnames(Sigma_NW))) {
    Sigma_NW <- Sigma_NW[var_order, var_order, drop = FALSE]
  } else {
    dimnames(Sigma_NW) <- list(var_order, var_order)
  }

  bm_ess   <- matrix_to_bmatrix(Sigma_ess,   digits = digits)
  bm_naive <- matrix_to_bmatrix(Sigma_naive, digits = digits)
  bm_NW    <- matrix_to_bmatrix(Sigma_NW,    digits = digits)

  cat("\\begin{table}[h]\n")
  cat("\\centering\n")
  cat("\\caption{Multivariate Shewhart individuals covariance matrices. Variable order is CO$_2$, humidity, soil temperature. Estimation based on ESS variance and covariance plots, the naive estimator corresponding to $T=1$, and the Newey--West estimator.}\n")
  cat("\\label{tab:Sigma-greenhouse}\n")
  cat("\\renewcommand{\\arraystretch}{1.15}\n")
  cat("\\begin{tabular}{ccc}\n")
  cat("\\toprule\n")
  cat("ESS & Naive & Newey--West \\\\\n")
  cat("\\midrule\n")
  cat(bm_ess, " &\n", bm_naive, " &\n", bm_NW, "\\\\\n", sep = "")
  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
}

print_greenhouse_Sigma_table(
  res_T,
  Sigma_NW   = Sigma_NW,
  var_order  = c("co2_ppm", "humidity_pct", "soil_temp_F"),
  method     = "roll",
  digits     = 1
)


