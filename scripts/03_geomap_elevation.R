suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(terra)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(scales)
  library(magick)
})

# -----------------------------
# Paths
# -----------------------------
geoloc_dir <- "data/geoloc"
out_dir <- "data/Elevation"
fig1_dir <- "Figs/Fig1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "extra"), recursive = TRUE, showWarnings = FALSE)
dir.create(fig1_dir, recursive = TRUE, showWarnings = FALSE)

elevation_tif <- "/Users/nirwantandukar/Documents/Research/data/P_prediction/predictors/ELEVATION.tif"
if (!file.exists(elevation_tif)) {
  stop("Elevation tif not found: ", elevation_tif)
}
elevation <- terra::rast(elevation_tif)

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

mix_colors <- function(col_a, col_b, w = 0.5) {
  a <- grDevices::col2rgb(col_a) / 255
  b <- grDevices::col2rgb(col_b) / 255
  m <- (1 - w) * a + w * b
  grDevices::rgb(m[1], m[2], m[3])
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
  out$elevation_m <- terra::extract(rast_obj, terra::vect(pts_sf))[, 2]
  out
}

# -----------------------------
# Load geoloc
# -----------------------------
at <- clean_names(read.csv(file.path(geoloc_dir, "at_1001_geoloc.csv"), check.names = FALSE))
barley <- clean_names(read.csv(file.path(geoloc_dir, "barley_geoloc.csv"), check.names = FALSE))
maize <- clean_names(read.csv(file.path(geoloc_dir, "maize_SEEDS_geoloc.csv"), check.names = FALSE))
rice <- clean_names(read.csv(file.path(geoloc_dir, "rice_3000_geoloc.csv"), check.names = FALSE))
sorghum <- clean_names(read.csv(file.path(geoloc_dir, "sorghum_Lasky_geoloc.csv"), check.names = FALSE))
pearl_path <- file.path(geoloc_dir, "extra", "pearl_millet_geoloc.csv")
has_pearl_millet <- file.exists(pearl_path)
if (has_pearl_millet) {
  pearl_millet <- clean_names(read.csv(pearl_path, check.names = FALSE))
}

at_e <- extract_values(at, id_col = "Acession_ID", lon_col = "Long", lat_col = "Lat", crop_name = "arabidopsis", rast_obj = elevation)
barley_e <- extract_values(barley, id_col = "Taxa", lon_col = "lon", lat_col = "lat", crop_name = "barley", rast_obj = elevation)
maize_e <- extract_values(maize, id_col = "Taxa", lon_col = "lon", lat_col = "lat", crop_name = "maize", rast_obj = elevation)
rice_e <- extract_values(rice, id_col = "Taxa", lon_col = "lon", lat_col = "lat", crop_name = "rice", rast_obj = elevation)
sorghum_e <- extract_values(sorghum, id_col = "Taxa", lon_col = "lon", lat_col = "lat", crop_name = "sorghum", rast_obj = elevation)

all_e <- bind_rows(maize_e, sorghum_e, rice_e, at_e, barley_e)
if (has_pearl_millet) {
  pearl_millet_e <- extract_values(pearl_millet, id_col = "Taxa", lon_col = "lon", lat_col = "lat", crop_name = "pearl_millet", rast_obj = elevation)
  all_e <- bind_rows(all_e, pearl_millet_e)
}
points <- all_e %>% filter(!is.na(elevation_m))
map_crop_order <- c("arabidopsis", "barley", "rice", "maize", "sorghum")
points <- points %>% filter(crop %in% map_crop_order)

# -----------------------------
# Save tables
# -----------------------------
write.csv(at_e, file.path(out_dir, "arabidopsis_elevation_values.csv"), row.names = FALSE)
write.csv(barley_e, file.path(out_dir, "barley_elevation_values.csv"), row.names = FALSE)
write.csv(maize_e, file.path(out_dir, "maize_elevation_values.csv"), row.names = FALSE)
write.csv(rice_e, file.path(out_dir, "rice_elevation_values.csv"), row.names = FALSE)
write.csv(sorghum_e, file.path(out_dir, "sorghum_elevation_values.csv"), row.names = FALSE)
if (has_pearl_millet) {
  write.csv(pearl_millet_e, file.path(out_dir, "extra", "pearl_millet_elevation_values.csv"), row.names = FALSE)
}
write.csv(all_e, file.path(out_dir, "all_crops_elevation_values.csv"), row.names = FALSE)

# -----------------------------
# Style
# -----------------------------
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
  maize = "#E69F00",
  sorghum = "#D55E00",
  rice = "#CC79A7",
  arabidopsis = "#009E73",
  barley = "#0072B2",
  pearl_millet = "#F0E442"
)

elev_min <- min(points$elevation_m, na.rm = TRUE)
elev_max <- max(points$elevation_m, na.rm = TRUE)
elev_range <- if (is.finite(elev_min) && is.finite(elev_max) && elev_max > elev_min) c(elev_min, elev_max) else c(0, 1)

points$col <- mapply(
  function(cr, val) {
    base_col <- crop_palettes[[cr]]
    # black-white style elevation shading while keeping crop hue cue
    low_col <- mix_colors(base_col, "#FFFFFF", w = 0.90)
    high_col <- mix_colors(base_col, "#000000", w = 0.55)
    col_numeric(c(low_col, high_col), domain = elev_range)(val)
  },
  points$crop,
  points$elevation_m
)

# -----------------------------
# Map
# -----------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

p <- ggplot() +
  geom_sf(data = world, fill = "white", color = "black", linewidth = 0.2) +
  geom_point(data = points, aes(x = Long, y = Lat), color = points$col, size = 0.6, alpha = 0.85)

crop_present <- map_crop_order[map_crop_order %in% unique(points$crop)]

legend_df <- data.frame(
  crop = crop_present,
  lon = -170,
  lat = seq(72, by = -7, length.out = length(crop_present)),
  col = unname(unlist(crop_palettes[crop_present]))
)

p <- p +
  annotate("rect", xmin = -178, xmax = -126, ymin = 20, ymax = 80, fill = alpha("white", 0.9), color = "black", linewidth = 0.25) +
  annotate("text", x = -175, y = 77, label = "Crop", hjust = 0, size = 3.4, fontface = "bold") +
  geom_point(data = legend_df, aes(x = lon, y = lat), color = legend_df$col, size = 3) +
  geom_text(data = legend_df, aes(x = lon, y = lat, label = crop), hjust = 0, nudge_x = 6, size = 3.5)

# Grayscale elevation legend (high darker, low lighter)
grad_vals <- seq(0, 1, length.out = 60)
grad_df <- data.frame(
  lon = -170,
  lat = seq(12, -26, length.out = length(grad_vals)),
  col = gray(seq(0.0, 1.0, length.out = length(grad_vals)))
)

p <- p +
  annotate("rect", xmin = -178, xmax = -126, ymin = -30, ymax = 18, fill = alpha("white", 0.9), color = "black", linewidth = 0.25) +
  annotate("text", x = -175, y = 15, label = "Elevation (m)", hjust = 0, size = 3.4, fontface = "bold") +
  geom_point(data = grad_df, aes(x = lon, y = lat), color = grad_df$col, size = 3.6) +
  annotate("text", x = -164.5, y = 12, label = paste0("high (", round(elev_max, 0), ")"), hjust = 0, size = 3.0) +
  annotate("text", x = -164.5, y = -26, label = paste0("low (", round(elev_min, 0), ")"), hjust = 0, size = 3.0) +
  labs(x = NULL, y = NULL) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE, datum = NA) +
  plot_theme +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

map_out_data <- file.path(out_dir, "geomap_elevation.png")
map_out_fig1 <- file.path(fig1_dir, "geomap_elevation.png")

ggsave(map_out_data, p, width = 12, height = 6, dpi = 200, device = ragg::agg_png, bg = "white")
ggsave(map_out_fig1, p, width = 12, height = 6, dpi = 200, device = ragg::agg_png, bg = "white")

# -----------------------------
# Fig1 assembly (timeline on top, map below)
# -----------------------------
timeline_png <- file.path(fig1_dir, "timeline_crops.png")
fig1_png <- file.path(fig1_dir, "Fig1_timeline_map.png")

if (file.exists(timeline_png) && file.exists(map_out_fig1)) {
  img_timeline <- magick::image_read(timeline_png)
  img_map <- magick::image_read(map_out_fig1)

  target_width <- max(
    magick::image_info(img_timeline)$width,
    magick::image_info(img_map)$width
  )

  img_timeline <- magick::image_resize(img_timeline, paste0(target_width, "x"))
  img_map <- magick::image_resize(img_map, paste0(target_width, "x"))

  info_t <- magick::image_info(img_timeline)
  timeline_h <- as.integer(info_t$height[[1]])

  fig1 <- magick::image_append(c(img_timeline, img_map), stack = TRUE) %>%
    magick::image_annotate("A", location = "+40+35", size = 120, weight = 700, color = "black") %>%
    magick::image_annotate("B", location = paste0("+40+", timeline_h + 35), size = 120, weight = 700, color = "black")
  magick::image_write(fig1, path = fig1_png, format = "png")
  message("Saved Fig1: ", fig1_png)
} else {
  warning("Could not build Fig1. Missing file: ", timeline_png)
}

message("Saved elevation map and tables to: ", out_dir)
