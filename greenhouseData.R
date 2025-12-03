rm(list=ls())
require(tidyverse)
require(lubridate)
require(patchwork)
require(xtable)
require(zoo)

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

#write data for repository
write_csv(phase1, file="GreenhouseEnvironmentTimeSeries.csv")
