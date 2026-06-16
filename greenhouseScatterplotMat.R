plot_greenhouse_scatter <- function(data,
                                        col = gh_cols("moody"),
                                        grid_col = gh_cols("grey_light"),
                                        T = 1,
                                        alphaT=.35) {

  vars <- c("co2_ppm", "humidity_pct", "soil_temp_F")
  stopifnot(all(vars %in% names(data)))

  # Order by time if available (so "every Tth" is meaningful)
  if ("date" %in% names(data)) {
    data <- data[order(data$date), , drop = FALSE]
  }

  # Keep every Tth observation
  if (!is.numeric(T) || length(T) != 1 || is.na(T) || T < 1) stop("T must be a positive number >= 1")
  T <- as.integer(T)
  idx <- seq.int(1, nrow(data), by = T)
  data_s <- data[idx, , drop = FALSE]

  X <- data_s[vars]

  pretty <- list(
    co2_ppm      = expression(CO[2]~"(ppm)"),
    humidity_pct = "Humidity (%)",
    soil_temp_F  = expression("Soil Temp ("*degree*"F)")
  )

  pt_col <- adjustcolor(col, alpha.f = alphaT)

  panel_scatter <- function(x, y, xlab, ylab) {
    plot(x, y,
         xlab = xlab, ylab = ylab,
         pch = 16, cex = 0.7, col = pt_col,
         axes = TRUE, frame.plot = TRUE)

    abline(h = axTicks(2), v = axTicks(1),
           col = grid_col, lty = "dashed", lwd = 0.7)

    points(x, y, pch = 16, cex = 0.7, col = pt_col)

    box(col = "grey30")
  }

  op <- par(
    bg = "white",
    fg = "grey20",
    col.axis = "grey20",
    col.lab = "grey20",
    family = "sans",
    las = 1,
    xaxs = "i", yaxs = "i",
    mar = c(4, 4, 2, 1) + 0.1,
    mgp = c(2.2, 0.8, 0),
    tcl = -0.25
  )
  on.exit(par(op), add = TRUE)

  # [1] A-B   [2] blank
  # [3] A-C   [4] B-C
  layout(matrix(c(1, 2,
                  3, 4), 2, 2, byrow = TRUE))

  panel_scatter(X$co2_ppm,      X$humidity_pct, pretty$co2_ppm,      pretty$humidity_pct)
  plot.new()  # blank upper-right
  panel_scatter(X$co2_ppm,      X$soil_temp_F,  pretty$co2_ppm,      pretty$soil_temp_F)
  panel_scatter(X$humidity_pct, X$soil_temp_F,  pretty$humidity_pct, pretty$soil_temp_F)

  invisible(NULL)
}

plot_greenhouse_scatter_24h <- function(data,
                                        start_time = NULL,   # e.g. "2025-08-12 12:01:00"
                                        day = NULL,          # e.g. as.Date("2025-08-12") or "2025-08-12"
                                        tz = attr(data$date, "tzone"),
                                        hours = 24,
                                        ...                  # passed to plot_greenhouse_scatter (col, T, alphaT, etc.)
                                        ) {

  stopifnot("date" %in% names(data))
  stopifnot(inherits(data$date, c("POSIXct", "POSIXt")))

  # choose window start
  if (!is.null(day)) {
    day <- as.Date(day)
    start <- as.POSIXct(day, tz = tz)
  } else if (!is.null(start_time)) {
    start <- as.POSIXct(start_time, tz = tz)
  } else {
    start <- min(data$date, na.rm = TRUE)
  }
  end <- start + hours * 3600

  data_24 <- data[data$date >= start & data$date < end, , drop = FALSE]
  if (nrow(data_24) == 0) stop("No rows found in the selected 24-hour window.")

  plot_greenhouse_scatter(data_24, ...)
  invisible(list(start = start, end = end, n = nrow(data_24)))
}

# Example: keep every 50th observation
pdf("greenhouse_scatterT50.pdf", width = 8, height = 8)
plot_greenhouse_scatter(phase1, col = gh_cols("moody"), T = 50)
dev.off()

pdf("greenhouse_scatterT1.pdf", width = 8, height = 8)
plot_greenhouse_scatter(phase1, col = gh_cols("moody"), T = 1, alphaT=.002)
dev.off()

pdf("greenhouse_scatterT1notrans.pdf", width = 8, height = 8)
plot_greenhouse_scatter(phase1, col = gh_cols("moody"), T = 1, alphaT=1)
dev.off()

pdf("greenhouse_scatter_2025-09-06.pdf", width = 8, height = 8)
plot_greenhouse_scatter_24h(
  phase1,
  day = "2025-09-06",
  col = gh_cols("moody"),
  T = 1,
  alphaT = .85
)
dev.off()
