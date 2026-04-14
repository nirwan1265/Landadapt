suppressPackageStartupMessages({
  library(dplyr)
  library(terra)
  library(geodata)
})

# ------------------------------------------------------------
# Paths and settings
# ------------------------------------------------------------
geoloc_dir <- "data/geoloc"
out_dir <- "data/WorldClim"
cache_dir <- file.path(out_dir, "cache")

# WorldClim resolution in arc-minutes (10 is lighter/faster)
worldclim_res <- 10

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
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
  hit <- names(df)[grepl(col_name, names(df), fixed = TRUE)]
  if (length(hit) > 0) return(hit[1])
  stop("Required column not found: ", col_name)
}

read_crop_points <- function(file_name, id_col, lon_col, lat_col, crop_name) {
  f <- file.path(geoloc_dir, file_name)
  if (!file.exists(f)) {
    message("Skipping missing file: ", f)
    return(NULL)
  }

  df <- clean_names(read.csv(f, check.names = FALSE))
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
    filter(!is.na(Long), !is.na(Lat)) %>%
    filter(Long >= -180, Long <= 180, Lat >= -90, Lat <= 90)

  out
}

# ------------------------------------------------------------
# Load crop geolocations (all current crops)
# ------------------------------------------------------------
crop_cfg <- list(
  list(file = "at_1001_geoloc.csv", id = "Acession_ID", lon = "Long", lat = "Lat", crop = "arabidopsis"),
  list(file = "barley_geoloc.csv", id = "Taxa", lon = "lon", lat = "lat", crop = "barley"),
  list(file = "maize_SEEDS_geoloc.csv", id = "Taxa", lon = "lon", lat = "lat", crop = "maize"),
  list(file = "rice_3000_geoloc.csv", id = "Taxa", lon = "lon", lat = "lat", crop = "rice"),
  list(file = "sorghum_Lasky_geoloc.csv", id = "Taxa", lon = "lon", lat = "lat", crop = "sorghum")
)

points_list <- lapply(crop_cfg, function(cfg) {
  read_crop_points(cfg$file, cfg$id, cfg$lon, cfg$lat, cfg$crop)
})

points <- bind_rows(points_list)
if (nrow(points) == 0) {
  stop("No geolocation points found. Check files in data/geoloc.")
}

# ------------------------------------------------------------
# Get WorldClim BIO01..BIO19 and extract at points
# ------------------------------------------------------------
message("Downloading/reading WorldClim BIO variables...")
bio <- geodata::worldclim_global(var = "bio", res = worldclim_res, path = cache_dir)

if (nlyr(bio) < 19) {
  stop("WorldClim bio raster has fewer than 19 layers.")
}

bio <- bio[[1:19]]
names(bio) <- sprintf("BIO%02d", 1:19)

bio_vals <- terra::extract(bio, points[, c("Long", "Lat")], ID = FALSE)
all_df <- bind_cols(points, bio_vals)

# ------------------------------------------------------------
# Crop-specific PCA (PC1..PC3 from BIO01..BIO19)
# ------------------------------------------------------------
bio_cols <- sprintf("BIO%02d", 1:19)
all_df$PC1 <- NA_real_
all_df$PC2 <- NA_real_
all_df$PC3 <- NA_real_

for (cr in unique(all_df$crop)) {
  idx <- which(all_df$crop == cr)
  sub <- all_df[idx, c("Lines", bio_cols)]
  keep <- complete.cases(sub[, bio_cols])

  if (sum(keep) < 3) {
    message("Skipping PCA for ", cr, " (fewer than 3 complete rows).")
    next
  }

  mat <- as.matrix(sub[keep, bio_cols])
  pca <- prcomp(mat, center = TRUE, scale. = TRUE)

  n_pc <- min(3, ncol(pca$x))
  pcs <- as.data.frame(pca$x[, seq_len(n_pc), drop = FALSE])
  names(pcs) <- paste0("PC", seq_len(n_pc))

  # Fill available PCs; keep missing PCs as NA
  keep_idx <- idx[keep]
  for (k in seq_len(n_pc)) {
    all_df[keep_idx, paste0("PC", k)] <- pcs[[k]]
  }

  # Save PCA loadings and variance explained per crop
  loadings <- as.data.frame(pca$rotation[, seq_len(n_pc), drop = FALSE])
  loadings$Variable <- rownames(loadings)
  loadings <- loadings %>% select(Variable, everything())
  write.csv(loadings, file.path(out_dir, paste0("worldclim_", cr, "_PCA_loadings.csv")), row.names = FALSE)

  var_exp <- data.frame(
    PC = paste0("PC", seq_along(pca$sdev)),
    VarianceExplained = (pca$sdev^2) / sum(pca$sdev^2)
  )
  write.csv(var_exp, file.path(out_dir, paste0("worldclim_", cr, "_PCA_variance.csv")), row.names = FALSE)
}

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
for (cr in unique(all_df$crop)) {
  out_file <- file.path(out_dir, paste0("worldclim_", cr, "_BIO_PCA.csv"))
  write.csv(all_df %>% filter(crop == cr), out_file, row.names = FALSE)
}

write.csv(all_df, file.path(out_dir, "worldclim_all_crops_BIO_PCA.csv"), row.names = FALSE)

message("Done. Saved WorldClim BIO and PCA files to: ", out_dir)
