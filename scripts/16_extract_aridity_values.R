suppressPackageStartupMessages({
  library(dplyr)
  library(terra)
  library(sf)
})

# -----------------------------
# Paths
# -----------------------------
geoloc_dir <- "data/geoloc"
out_dir <- "data/Aridity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "extra"), recursive = TRUE, showWarnings = FALSE)

ai_tif <- "/Users/nirwantandukar/Documents/Research/data/P_prediction/Global-AI_ET0_v3_annual/ai_v3_yr.tif"
if (!file.exists(ai_tif)) stop("Aridity tif not found: ", ai_tif)

ai_rast <- terra::rast(ai_tif)

# Dataset uses 65535 as nodata sentinel.
terra::NAflag(ai_rast) <- 65535

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

extract_values <- function(df, id_col, lon_col, lat_col, crop_name, rast_obj) {
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

  # Returns scaled integer values in most AI v3 distributions.
  out$aridity_index_raw <- terra::extract(rast_obj, terra::vect(pts_sf))[, 2]
  out$aridity_index <- ifelse(
    is.na(out$aridity_index_raw),
    NA_real_,
    as.numeric(out$aridity_index_raw) / 10000
  )

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
# Run extraction
# -----------------------------
all_vals <- list()
for (x in cfg) {
  df <- read_geo(x$file)
  out <- extract_values(df, x$id, x$lon, x$lat, x$name, ai_rast)

  if (x$name == "pearl_millet") {
    out_file <- file.path(out_dir, "extra", "pearl_millet_aridity_values.csv")
  } else {
    out_file <- file.path(out_dir, paste0(x$name, "_aridity_values.csv"))
  }
  write.csv(out, out_file, row.names = FALSE)
  all_vals[[x$name]] <- out
}

all_df <- bind_rows(all_vals)
write.csv(all_df, file.path(out_dir, "all_crops_aridity_values.csv"), row.names = FALSE)

message("Saved aridity values to: ", out_dir)
message("Columns: Lines, Long, Lat, crop, aridity_index_raw, aridity_index")
