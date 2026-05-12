suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(patchwork)
  library(scales)
  library(png)
  library(grid)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
data_dir <- "data/WorldClim"
res_dir <- "results"
fig_dir <- "Figs/Supplementary"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Panel A: species-level scree mean +/- SE
# ------------------------------------------------------------
var_files <- list.files(
  data_dir,
  pattern = "^worldclim_(arabidopsis|barley|rice|maize|sorghum)_PCA_variance\\.csv$",
  full.names = TRUE
)
if (length(var_files) == 0) {
  stop("No species PCA variance files found in ", data_dir)
}

var_df <- bind_rows(lapply(var_files, function(path) {
  sp <- basename(path) |>
    str_remove("^worldclim_") |>
    str_remove("_PCA_variance\\.csv$")

  read_csv(path, show_col_types = FALSE) |>
    mutate(species = sp)
}))

scree_species <- var_df |>
  mutate(
    pc_num = as.integer(str_remove(PC, "^PC")),
    var_pct = 100 * VarianceExplained
  ) |>
  filter(pc_num <= 10)

scree_mean <- scree_species |>
  group_by(pc_num) |>
  summarise(
    mean_var = mean(var_pct, na.rm = TRUE),
    se_var = sd(var_pct, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(PC = factor(paste0("PC", pc_num), levels = paste0("PC", 1:10)))

pA <- ggplot(scree_mean, aes(x = PC, y = mean_var, group = 1)) +
  geom_col(fill = "grey35", width = 0.72) +
  geom_errorbar(
    aes(ymin = pmax(mean_var - se_var, 0), ymax = mean_var + se_var),
    width = 0.2,
    linewidth = 0.35
  ) +
  geom_line(color = "black", linewidth = 0.4) +
  geom_point(color = "black", size = 1.7) +
  labs(
    title = "Shared PCA Scree Across Species",
    x = NULL,
    y = "Mean variance explained (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

# ------------------------------------------------------------
# Panel B: pooled PCA separation + loading arrows
# ------------------------------------------------------------
bio_files <- list.files(
  data_dir,
  pattern = "^worldclim_(arabidopsis|barley|rice|maize|sorghum)_BIO_PCA\\.csv$",
  full.names = TRUE
)
if (length(bio_files) == 0) {
  stop("No species BIO PCA source files found in ", data_dir)
}

bio_cols <- sprintf("BIO%02d", 1:19)

pooled <- bind_rows(lapply(bio_files, function(path) {
  sp <- basename(path) |>
    str_remove("^worldclim_") |>
    str_remove("_BIO_PCA\\.csv$")
  dat <- read_csv(path, show_col_types = FALSE)
  if (!all(bio_cols %in% names(dat))) {
    stop("Missing BIO01..BIO19 in: ", path)
  }
  dat |>
    mutate(species = sp) |>
    select(species, all_of(bio_cols))
}))

pooled_clean <- pooled |>
  filter(if_all(all_of(bio_cols), ~ !is.na(.x)))

if (nrow(pooled_clean) < 20) {
  stop("Too few complete rows for pooled PCA")
}

pca <- prcomp(as.matrix(pooled_clean[, bio_cols]), center = TRUE, scale. = TRUE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as_tibble(pca$x) |>
  bind_cols(pooled_clean |> select(species)) |>
  mutate(
    species = factor(species, levels = c("arabidopsis", "barley", "rice", "maize", "sorghum"))
  )

loadings <- as_tibble(pca$rotation, rownames = "Variable") |>
  mutate(vec_norm = sqrt(PC1^2 + PC2^2 + PC3^2)) |>
  arrange(desc(vec_norm)) |>
  slice_head(n = 10)

species_cols <- c(
  arabidopsis = "#1B9E77",
  barley = "#7570B3",
  rice = "#E7298A",
  maize = "#D95F02",
  sorghum = "#66A61E"
)

make_pair_plot <- function(xpc, ypc) {
  xq <- quantile(abs(scores[[xpc]]), 0.95, na.rm = TRUE)
  yq <- quantile(abs(scores[[ypc]]), 0.95, na.rm = TRUE)
  arrow_mult <- 0.72 * min(
    xq / max(abs(loadings[[xpc]]), na.rm = TRUE),
    yq / max(abs(loadings[[ypc]]), na.rm = TRUE)
  )

  arrow_df <- loadings |>
    transmute(
      Variable,
      x = 0,
      y = 0,
      xend = .data[[xpc]] * arrow_mult,
      yend = .data[[ypc]] * arrow_mult
    )

  xnum <- as.integer(str_remove(xpc, "^PC"))
  ynum <- as.integer(str_remove(ypc, "^PC"))

  ggplot(scores, aes(x = .data[[xpc]], y = .data[[ypc]], color = species)) +
    geom_point(alpha = 0.6, size = 1.2, stroke = 0) +
    stat_ellipse(level = 0.68, linewidth = 0.35, alpha = 0.9, show.legend = FALSE) +
    geom_segment(
      data = arrow_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.010, "npc")),
      color = "grey20",
      linewidth = 0.3
    ) +
    geom_text(
      data = arrow_df,
      aes(x = xend, y = yend, label = Variable),
      inherit.aes = FALSE,
      color = "grey10",
      size = 2.2,
      vjust = -0.2
    ) +
    scale_color_manual(values = species_cols) +
    labs(
      x = paste0(xpc, " (", sprintf("%.1f", 100 * var_expl[xnum]), "%)"),
      y = paste0(ypc, " (", sprintf("%.1f", 100 * var_expl[ynum]), "%)"),
      color = "Species"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
}

pB12 <- make_pair_plot("PC1", "PC2")
pB13 <- make_pair_plot("PC1", "PC3")
pB23 <- make_pair_plot("PC2", "PC3")

pB <- (pB12 | pB13 | pB23) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pB <- pB + plot_annotation(title = "Shared PCA Separation and Loading Vectors (PC1, PC2, PC3)")

# ------------------------------------------------------------
# Panel C + D: pre-rendered heatmaps
# ------------------------------------------------------------
c_path <- file.path(fig_dir, "WorldClim_PCA_loadings_species_rows_heatmap.png")
d_path <- file.path(fig_dir, "WorldClim_PCA_meanAbs_loadings_heatmap.png")
if (!file.exists(c_path) || !file.exists(d_path)) {
  stop("Missing heatmap input(s): ", c_path, " and/or ", d_path)
}

c_img <- rasterGrob(readPNG(c_path), interpolate = TRUE)
d_img <- rasterGrob(readPNG(d_path), interpolate = TRUE)

pC <- wrap_elements(full = c_img) + theme_void()
pD <- wrap_elements(full = d_img) + theme_void()

# ------------------------------------------------------------
# Compose and save
# ------------------------------------------------------------
fig <- ((pA | pB) + plot_layout(widths = c(1, 2))) / (pC | pD) +
  plot_annotation(tag_levels = "A")

out_png <- file.path(fig_dir, "SuppFig_WorldClim_PCA.png")
out_pdf <- file.path(fig_dir, "SuppFig_WorldClim_PCA.pdf")

ggsave(out_png, fig, width = 18, height = 14, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(out_pdf, fig, width = 18, height = 14, bg = "white", device = cairo_pdf)

message("Saved: ", out_png)
message("Saved: ", out_pdf)
