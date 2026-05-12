suppressPackageStartupMessages({
  library(dplyr)
  library(terra)
  library(sf)
})

# -----------------------------
# Paths
# -----------------------------
geoloc_dir <- "data/geoloc"
predictor_dir <- "/Users/nirwantandukar/Documents/Research/data/P_prediction/predictors"

pick_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

matemp_tif <- pick_first_existing(c(
  file.path(predictor_dir, "MATEMP.tif"),
  file.path(predictor_dir, "wc2.1_30s_bio_1.tif"),
  file.path(predictor_dir, "root biomass carbon", "wc2.1_30s_bio_1.tif")
))
maprecip_tif <- pick_first_existing(c(
  file.path(predictor_dir, "MAPRECIP.tif")
))

layers <- list(
  list(
    key = "MATEMP",
    tif = matemp_tif,
    out_dir = "data/MATEMP",
    value_col = "matemp_value"
  ),
  list(
    key = "MAPRECIP",
    tif = maprecip_tif,
    out_dir = "data/MAPRECIP",
    value_col = "maprecip_value"
  )
)

for (x in layers) {
  if (is.na(x$tif)) stop("Missing tif for layer: ", x$key)
  dir.create(x$out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(x$out_dir, "extra"), recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------
# Helpers
# -----------------------------
clean_names <- function(df) {
  nm <- names(df)
  nm <- sub("^\\ufeff", "", nm, useBytes = TRUE)
  bad <- nm == "" | is.na(nm)
  df <- df[, !bad, drop = FALSE]
  names(df) <- nm[!bad]
  df
}

resolve_col <- function(df, col_name) {
  if (col_name %in% names(df)) return(col_name)
  hits <- names(df)[grepl(col_name, names(df), fixed = TRUE)]
  if (length(hits) > 0) return(hits[1])
  stop("Required column not found: ", col_name)
}

extract_values <- function(df, id_col, lon_col, lat_col, crop_name, rast_obj, value_col) {
  id_col <- resolve_col(df, id_col)
  lon_col <- resolve_col(df, lon_col)
  lat_col <- resolve_col(df, lat_col)

  out <- df %>%
    transmute(
      Lines = as.character(.data[[id_col]]),
      Long = as.numeric(.data[[lon_col]]),
      Lat = as.numeric(.data[[lat_col]]),
      crop = crop_name
    ) %>%
    filter(!is.na(Long), !is.na(Lat))

  pts_sf <- st_as_sf(out, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)
  pts_sf <- st_transform(pts_sf, crs(rast_obj))
  out[[value_col]] <- terra::extract(rast_obj, terra::vect(pts_sf))[, 2]
  out
}

read_geo <- function(path) clean_names(read.csv(path, check.names = FALSE))

# -----------------------------
# Geolocation config
# -----------------------------
cfg <- list(
  list(name = "arabidopsis", file = file.path(geoloc_dir, "at_1001_geoloc.csv"), id = "Acession_ID", lon = "Long", lat = "Lat"),
  list(name = "barley", file = file.path(geoloc_dir, "barley_geoloc.csv"), id = "Taxa", lon = "lon", lat = "lat"),
  list(name = "maize", file = file.path(geoloc_dir, "maize_SEEDS_geoloc.csv"), id = "Taxa", lon = "lon", lat = "lat"),
  list(name = "rice", file = file.path(geoloc_dir, "rice_3000_geoloc.csv"), id = "Taxa", lon = "lon", lat = "lat"),
  list(name = "sorghum", file = file.path(geoloc_dir, "sorghum_Lasky_geoloc.csv"), id = "Taxa", lon = "lon", lat = "lat")
)

pearl_path <- file.path(geoloc_dir, "extra", "pearl_millet_geoloc.csv")
if (file.exists(pearl_path)) {
  cfg <- append(cfg, list(list(name = "pearl_millet", file = pearl_path, id = "Taxa", lon = "lon", lat = "lat")))
}

# -----------------------------
# Run extraction for all layers
# -----------------------------
for (layer in layers) {
  r <- terra::rast(layer$tif)
  all_vals <- list()

  for (x in cfg) {
    df <- read_geo(x$file)
    out <- extract_values(df, x$id, x$lon, x$lat, x$name, r, layer$value_col)

    if (x$name == "pearl_millet") {
      out_file <- file.path(layer$out_dir, "extra", paste0("pearl_millet_", layer$key, "_values.csv"))
    } else {
      out_file <- file.path(layer$out_dir, paste0(x$name, "_", layer$key, "_values.csv"))
    }
    write.csv(out, out_file, row.names = FALSE)
    all_vals[[x$name]] <- out
  }

  all_df <- bind_rows(all_vals)
  write.csv(all_df, file.path(layer$out_dir, paste0("all_crops_", layer$key, "_values.csv")), row.names = FALSE)
  message("Saved layer values: ", layer$key, " (", layer$tif, ") -> ", layer$out_dir)
}

