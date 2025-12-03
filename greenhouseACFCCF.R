# ACF and CCF -------------------------------------------------------------

# settings
lag_max <- 20000L
period  <- 1440L  # 24 hours in minutes

# labels as plotmath strings (for label_parsed)
acf_labels <- c(
  co2_ppm      = "CO[2]",
  humidity_pct = "Humidity",
  soil_temp_F  = "Soil~Temp"
)

ccf_pairs <- list(
  "co2_ppm|humidity_pct"      = c("co2_ppm",      "humidity_pct"),
  "humidity_pct|soil_temp_F"  = c("humidity_pct", "soil_temp_F"),
  "soil_temp_F|co2_ppm"       = c("soil_temp_F",  "co2_ppm")
)

ccf_labels <- c(
  "co2_ppm|humidity_pct"      = "CO[2]~'vs'~Humidity",
  "humidity_pct|soil_temp_F"  = "Humidity~'vs'~'Soil Temp'",
  "soil_temp_F|co2_ppm"       = "'Soil Temp'~'vs'~CO[2]"
)

# Maintain consistent row order
acf_levels <- unname(acf_labels)
ccf_levels <- unname(ccf_labels[names(ccf_pairs)])

# Helpers

acf_df <- function(x, lag_max) {
  a <- acf(x, lag.max = lag_max, plot = FALSE, na.action = na.pass)
  tibble(
    lag  = as.integer(a$lag),
    corr = as.numeric(a$acf)
  ) |>
    filter(lag >= 0)   # ACF: nonnegative lags only
}

ccf_df <- function(x, y, lag_max) {
  c <- ccf(x, y, lag.max = lag_max, plot = FALSE, na.action = na.pass)
  tibble(
    lag  = as.integer(c$lag),       # CCF: keep ± lags
    corr = as.numeric(c$acf)
  )
}

# Build tidy ACF for the left column 

acfs <- imap_dfr(acf_labels, ~ {
  var_name <- .y     # name in acf_labels
  label    <- .x     # plotmath string
  out <- acf_df(phase1[[var_name]], lag_max)
  out$panel_label <- label
  out
}) |>
  mutate(panel_label = factor(panel_label, levels = acf_levels))

vlines_acf <- tibble(lag = seq(0, lag_max, by = period))

p_acf <- ggplot(acfs, aes(lag, corr)) +
  geom_segment(aes(xend = lag, y = 0, yend = corr),
               linewidth = 0.35,
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
  facet_wrap(
    ~ panel_label,
    ncol = 1,
    scales = "fixed",
    strip.position = "top",
    labeller = label_parsed
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Lag", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    strip.text         = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Build tidy CCF for right column 

ccfs <- imap_dfr(ccf_pairs, ~ {
  var_x <- .x[1]
  var_y <- .x[2]
  key   <- .y
  out <- ccf_df(phase1[[var_x]], phase1[[var_y]], lag_max)
  out$panel_label <- ccf_labels[[key]]
  out
}) |>
  mutate(panel_label = factor(panel_label, levels = ccf_levels))

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
  facet_wrap(
    ~ panel_label,
    ncol = 1,
    scales = "fixed",
    strip.position = "top",
    labeller = label_parsed
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Lag", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    strip.text         = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Combine acfs and ccfs

greenhouseacfccfs <- p_acf | p_ccf

ggsave("acfccf.pdf", 
       greenhouseacfccfs,
       width=10,
       height=6,
       units="in")
