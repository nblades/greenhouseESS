plot_greenhouse_facets <- function(data, 
                                   col = gh_cols("moody"), 
                                   red_line = NULL) {
#red_line allows you to add a vertical line at a concerning date
  
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
