library(GGally)

pretty_labels <- c(
  co2_ppm      = "CO2 (ppm)",
  humidity_pct = "Humidity (%)",
  soil_temp_F  = "Soil Temp (°F)"
)

greenhouse_pairs <- ggpairs(
  phase1 |> select(co2_ppm, humidity_pct, soil_temp_F) |>
    rename("CO2 (ppm)"      = co2_ppm,
           "Humidity (%)"   = humidity_pct,
           "Soil Temp (°F)" = soil_temp_F),
  lower = list(continuous = wrap("points", color = gh_cols("moody"), alpha = 0.05, size = 0.3)),
  diag  = list(continuous = wrap("densityDiag", color = gh_cols("moody"), fill = gh_cols("moody"), alpha = 0.4)),
  upper = list(continuous = wrap("cor", size = 4, stars = FALSE))
) +
  theme_bw() +
  theme(
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.3),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "plain")
  )


ggsave("greenhouse_pairs.pdf", greenhouse_pairs)
