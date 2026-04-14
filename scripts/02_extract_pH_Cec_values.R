suppressPackageStartupMessages({
  library(dplyr)
  library(terra)
  library(sf)
})

# -----------------------------
# Paths
# -----------------------------
geoloc_dir <- "data/geoloc"
out_ph_dir <- "data/pH"
out_cec_dir <- "data/Cec"
dir.create(out_ph_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_cec_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_ph_dir, "extra"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_cec_dir, "extra"), recursive = TRUE, showWarnings = FALSE)

ph_tif <- "/Users/nirwantandukar/Documents/Research/data/P_prediction/predictors/phh2o_0-5cm_mean_1000.tif"
cec_tif <- "/Users/nirwantandukar/Documents/Research/data/P_prediction/predictors/cec_0-5cm_mean_1000.tif"

if (!file.exists(ph_tif)) stop("pH tif not found: ", ph_tif)
if (!file.exists(cec_tif)) stop("CEC tif not found: ", cec_tif)

ph_rast <- terra::rast(ph_tif)
cec_rast <- terra::rast(cec_tif)

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
# Load geolocation files
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

all_ph <- list()
all_cec <- list()

for (x in cfg) {
  df <- read_geo(x$file)
  ph_df <- extract_values(df, x$id, x$lon, x$lat, x$name, ph_rast, "ph_value")
  cec_df <- extract_values(df, x$id, x$lon, x$lat, x$name, cec_rast, "cec_value")

  if (x$name == "pearl_millet") {
    write.csv(ph_df, file.path(out_ph_dir, "extra", "pearl_millet_pH_values.csv"), row.names = FALSE)
    write.csv(cec_df, file.path(out_cec_dir, "extra", "pearl_millet_Cec_values.csv"), row.names = FALSE)
  } else {
    write.csv(ph_df, file.path(out_ph_dir, paste0(x$name, "_pH_values.csv")), row.names = FALSE)
    write.csv(cec_df, file.path(out_cec_dir, paste0(x$name, "_Cec_values.csv")), row.names = FALSE)
  }

  all_ph[[x$name]] <- ph_df
  all_cec[[x$name]] <- cec_df
}

all_ph_df <- bind_rows(all_ph)
all_cec_df <- bind_rows(all_cec)

write.csv(all_ph_df, file.path(out_ph_dir, "all_crops_pH_values.csv"), row.names = FALSE)
write.csv(all_cec_df, file.path(out_cec_dir, "all_crops_Cec_values.csv"), row.names = FALSE)

message("Saved pH values to: ", out_ph_dir)
message("Saved CEC values to: ", out_cec_dir)
