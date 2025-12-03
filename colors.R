# ---- Greenhouse unified colors ----
gh_palette <- list(
  moody      = "#36648B",  # primary line/stems (your greenhouse blue)
  ref        = "#63B8FF",  # reference (e.g., T=1 dashed)
  green      = "#009E73",  # 2nd categorical
  orange     = "#E69F00",  # 3rd categorical
  grey_light = "#D0D5DD",  # light grid/dash
  grey_mid   = "#9CA3AF",
  grey_dark  = "#6B7280"
)

gh_cols <- function(...) unname(unlist(gh_palette[c(...)]))

# consistent mapping for sample sizes n = {2000, 5000, 20000}
gh_scale_n <- function(name = "n") {
  scale_color_manual(
    values = c(`2000` = gh_cols("moody"),
               `5000` = gh_cols("orange"),
               `20000`= gh_cols("green")),
    limits = c("2000","5000","20000"),
    drop   = FALSE, name = name
  )
}

# a light, consistent base theme
theme_greenhouse <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background   = element_blank(),
      strip.placement    = "outside",
      strip.text.y.left  = element_text(angle = 0, hjust = 1, face = "plain")
    )
}
