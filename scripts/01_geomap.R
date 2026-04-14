library(dplyr)
library(ggplot2)
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)

# -----------------------------
# Paths
# -----------------------------
geoloc_dir <- "data/geoloc"
out_dir <- "data/TN"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nitrogen_tif <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/Phosphorus_prediction/data/phosphorus_prediction/predictors/nitrogen_0-5cm_mean_1000.tif"
if (!file.exists(nitrogen_tif)) {
  stop("Nitrogen tif not found: ", nitrogen_tif)
}
nitrogen <- terra::rast(nitrogen_tif)

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

extract_n_values <- function(df, id_col, lon_col, lat_col, crop_name, rast_obj) {
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

  # Reproject points from WGS84 to raster CRS before extraction.
  pts_sf <- st_as_sf(out, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)
  pts_sf <- st_transform(pts_sf, crs(rast_obj))

  out$n_value <- terra::extract(rast_obj, terra::vect(pts_sf))[, 2]
  out$log_n <- ifelse(out$n_value > 0, log10(out$n_value), NA_real_)
  out
}

# -----------------------------
# Load geolocation files (pearl millet optional)
# -----------------------------
at <- clean_names(read.csv(file.path(geoloc_dir, "at_1001_geoloc.csv"), check.names = FALSE))
barley <- clean_names(read.csv(file.path(geoloc_dir, "barley_geoloc.csv"), check.names = FALSE))
maize <- clean_names(read.csv(file.path(geoloc_dir, "maize_SEEDS_geoloc.csv"), check.names = FALSE))
rice <- clean_names(read.csv(file.path(geoloc_dir, "rice_3000_geoloc.csv"), check.names = FALSE))
sorghum <- clean_names(read.csv(file.path(geoloc_dir, "sorghum_Lasky_geoloc.csv"), check.names = FALSE))
pearl_millet_path <- file.path(geoloc_dir, "pearl_millet_geoloc.csv")
has_pearl_millet <- file.exists(pearl_millet_path)
if (has_pearl_millet) {
  pearl_millet <- clean_names(read.csv(pearl_millet_path, check.names = FALSE))
}

at_n <- extract_n_values(at, id_col = "Acession_ID", lon_col = "Long", lat_col = "Lat", crop_name = "arabidopsis", rast_obj = nitrogen)
barley_n <- extract_n_values(barley, id_col = "Taxa", lon_col = "lon", lat_col = "lat", crop_name = "barley", rast_obj = nitrogen)
maize_n <- extract_n_values(maize, id_col = "Taxa", lon_col = "lon", lat_col = "lat", crop_name = "maize", rast_obj = nitrogen)
rice_n <- extract_n_values(rice, id_col = "Taxa", lon_col = "lon", lat_col = "lat", crop_name = "rice", rast_obj = nitrogen)
sorghum_n <- extract_n_values(sorghum, id_col = "Taxa", lon_col = "lon", lat_col = "lat", crop_name = "sorghum", rast_obj = nitrogen)

all_n <- bind_rows(maize_n, sorghum_n, rice_n, at_n, barley_n)
if (has_pearl_millet) {
  pearl_millet_n <- extract_n_values(pearl_millet, id_col = "Taxa", lon_col = "lon", lat_col = "lat", crop_name = "pearl_millet", rast_obj = nitrogen)
  all_n <- bind_rows(all_n, pearl_millet_n)
}
points <- all_n
points <- points %>% filter(!is.na(log_n))

# Save TN tables from the same pipeline
write.csv(barley_n, file.path(out_dir, "barley_N_values.csv"), row.names = FALSE)
write.csv(maize_n, file.path(out_dir, "maize_N_values.csv"), row.names = FALSE)
if (has_pearl_millet) {
  write.csv(pearl_millet_n, file.path(out_dir, "pearl_millet_N_values.csv"), row.names = FALSE)
}
write.csv(sorghum_n, file.path(out_dir, "sorghum_N_values.csv"), row.names = FALSE)
write.csv(rice_n, file.path(out_dir, "rice_N_values.csv"), row.names = FALSE)
write.csv(at_n, file.path(out_dir, "arabidopsis_N_values.csv"), row.names = FALSE)
write.csv(all_n, file.path(out_dir, "all_crops_N_values.csv"), row.names = FALSE)

# -----------------------------
# Plot style and colors
# -----------------------------
split_by_crop <- FALSE

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  )

crop_palettes <- list(
  maize = c("#FFF6B0", "#D4A000"),
  sorghum = c("#F2D6B5", "#7A4A12"),
  rice = c("#FFD6D6", "#B30000"),
  arabidopsis = c("#CFEFC7", "#1B7A2A"),
  barley = c("#DCEBFF", "#2C6BB0"),
  pearl_millet = c("#F2D7FF", "#8E24AA")
)

log_min <- min(points$log_n, na.rm = TRUE)
log_max <- max(points$log_n, na.rm = TRUE)
log_range <- if (is.finite(log_min) && is.finite(log_max) && log_max > log_min) c(log_min, log_max) else c(0, 1)

points$col <- mapply(
  function(cr, val) {
    pal <- crop_palettes[[cr]]
    col_numeric(pal, domain = log_range)(val)
  },
  points$crop,
  points$log_n
)

# -----------------------------
# Map
# -----------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

p <- ggplot() +
  geom_sf(data = world, fill = "white", color = "black", linewidth = 0.2) +
  geom_point(data = points, aes(x = Long, y = Lat), color = points$col, size = 0.6, alpha = 0.85)

# Crop legend markers (South America ocean gap)
legend_df <- data.frame(
  crop = unique(points$crop),
  lon = -120,
  lat = seq(-6, by = -6, length.out = length(unique(points$crop))),
  col = sapply(
    unique(points$crop),
    function(cr) col_numeric(crop_palettes[[cr]], domain = log_range)(log_max)
  )
)

p <- p +
  geom_point(data = legend_df, aes(x = lon, y = lat), color = legend_df$col, size = 3) +
  geom_text(data = legend_df, aes(x = lon, y = lat, label = crop), hjust = 0, nudge_x = 6, size = 3.5)

# Vertical gradient legend for log10(N): white -> black
grad_vals <- seq(0, 1, length.out = 60)
grad_df <- data.frame(
  lon = -120,
  lat = seq(-46, -54, length.out = length(grad_vals)),
  col = gray(1 - grad_vals)
)

p <- p +
  geom_point(data = grad_df, aes(x = lon, y = lat), color = grad_df$col, size = 3.6) +
  annotate("text", x = -118.5, y = -46, label = "high", hjust = 0, size = 3.0) +
  annotate("text", x = -118.5, y = -54, label = "low", hjust = 0, size = 3.0) +
  annotate("text", x = -118.5, y = -50, label = "log10(N)", hjust = 0, size = 3.0) +
  coord_sf(expand = FALSE) +
  plot_theme +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

if (split_by_crop) {
  for (cr in unique(points$crop)) {
    p_cr <- ggplot() +
      geom_sf(data = world, fill = "white", color = "black", linewidth = 0.2) +
      geom_point(data = points %>% filter(crop == cr), aes(x = Long, y = Lat), color = points$col[points$crop == cr], size = 0.6, alpha = 0.85) +
      coord_sf(expand = FALSE) +
      plot_theme +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    ggsave(paste0("data/geomap_nitrogen_", cr, ".png"), p_cr, width = 12, height = 6, dpi = 200, device = ragg::agg_png, bg = "white")
  }
} else {
  ggsave("data/geomap_nitrogen.png", p, width = 12, height = 6, dpi = 200, device = ragg::agg_png, bg = "white")
}

# -----------------------------
# Histograms (log10(N), crop gradient colors)
# -----------------------------
for (cr in unique(points$crop)) {
  df_cr <- points %>% filter(crop == cr)
  if (nrow(df_cr) == 0) next
  pal <- crop_palettes[[cr]]

  h <- ggplot(df_cr, aes(x = log_n)) +
    geom_histogram(aes(fill = after_stat(x)), bins = 40, color = "white") +
    labs(
      title = paste("Nitrogen values:", cr),
      x = "log10(Nitrogen)",
      y = "Count"
    ) +
    scale_fill_gradient(low = pal[1], high = pal[2], name = "log10(N)") +
    plot_theme +
    theme(legend.position = "right")

  ggsave(paste0("data/nitrogen_hist_", cr, ".png"), h, width = 6, height = 4, dpi = 200, device = ragg::agg_png, bg = "white")
}

message("Saved map, histograms, and TN tables successfully.")
