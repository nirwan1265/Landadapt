suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
  library(png)
  library(grid)
  library(scatterplot3d)
})

data_dir <- "data/WorldClim"
fig_dir <- "Figs/Supplementary"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

species_order <- c("arabidopsis", "barley", "rice", "maize", "sorghum")
species_labels <- c(
  arabidopsis = "Arabidopsis",
  barley = "Barley",
  rice = "Rice",
  maize = "Maize",
  sorghum = "Sorghum"
)
species_cols <- c(
  arabidopsis = "#1B9E77",
  barley = "#7570B3",
  rice = "#E7298A",
  maize = "#D95F02",
  sorghum = "#66A61E"
)

# ------------------------------------------------------------
# SUPP FIG 1: A (scree mean +/- SE) + B (species-wise loadings heatmap)
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

scree_df <- var_df |>
  mutate(
    pc_num = as.integer(str_remove(PC, "^PC")),
    var_pct = 100 * VarianceExplained
  ) |>
  filter(pc_num <= 10) |>
  group_by(pc_num) |>
  summarise(
    mean_var = mean(var_pct, na.rm = TRUE),
    se_var = sd(var_pct, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(PC = factor(paste0("PC", pc_num), levels = paste0("PC", 1:10)))

pA <- ggplot(scree_df, aes(x = PC, y = mean_var, group = 1)) +
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

heatmap_path <- file.path(fig_dir, "WorldClim_PCA_loadings_species_rows_heatmap.png")
if (!file.exists(heatmap_path)) {
  stop("Missing file: ", heatmap_path)
}
pB <- wrap_elements(full = rasterGrob(readPNG(heatmap_path), interpolate = TRUE)) + theme_void()

suppfig1 <- (pA / pB) + plot_layout(heights = c(1, 1.8)) + plot_annotation(tag_levels = "A")

out1_png <- file.path(fig_dir, "SuppFig1_WorldClim_PCA.png")
out1_pdf <- file.path(fig_dir, "SuppFig1_WorldClim_PCA.pdf")
ggsave(out1_png, suppfig1, width = 14, height = 14, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(out1_pdf, suppfig1, width = 14, height = 14, bg = "white", device = cairo_pdf)

# ------------------------------------------------------------
# SUPP FIG 2: species-wise 3D PCA loading arrows (PC1, PC2, PC3)
# ------------------------------------------------------------
load_files <- list.files(
  data_dir,
  pattern = "^worldclim_(arabidopsis|barley|rice|maize|sorghum)_PCA_loadings\\.csv$",
  full.names = TRUE
)
if (length(load_files) == 0) {
  stop("No species PCA loading files found in ", data_dir)
}

top_n_arrows <- 12
load_by_species <- lapply(species_order, function(sp) {
  f <- load_files[grepl(paste0("worldclim_", sp, "_PCA_loadings\\.csv$"), load_files)]
  if (length(f) != 1) stop("Expected one PCA loading file for ", sp)
  dat <- read_csv(f, show_col_types = FALSE)
  req <- c("Variable", "PC1", "PC2", "PC3")
  if (!all(req %in% names(dat))) stop("Missing Variable/PC1/PC2/PC3 in file: ", f)
  dat |>
    transmute(
      Variable = as.character(Variable),
      PC1 = as.numeric(PC1),
      PC2 = as.numeric(PC2),
      PC3 = as.numeric(PC3)
    ) |>
    filter(if_all(c(PC1, PC2, PC3), ~ !is.na(.x))) |>
    mutate(vec_norm = sqrt(PC1^2 + PC2^2 + PC3^2)) |>
    arrange(desc(vec_norm)) |>
    slice_head(n = top_n_arrows)
})
names(load_by_species) <- species_order

all_load <- bind_rows(load_by_species)
lim_val <- max(abs(c(all_load$PC1, all_load$PC2, all_load$PC3)), na.rm = TRUE) * 1.08
x_lim <- c(-lim_val, lim_val)
y_lim <- c(-lim_val, lim_val)
z_lim <- c(-lim_val, lim_val)

draw_suppfig2 <- function(outfile, type = c("png", "pdf")) {
  type <- match.arg(type)
  if (type == "png") {
    ragg::agg_png(outfile, width = 3600, height = 2400, res = 300, background = "white")
  } else {
    cairo_pdf(outfile, width = 14, height = 9, bg = "white")
  }

  par(mfrow = c(2, 3), mar = c(3.2, 3.2, 2.5, 1.5), oma = c(0, 0, 1.5, 0))
  for (sp in species_order) {
    d <- load_by_species[[sp]]
    s3d <- scatterplot3d(
      x = d$PC1,
      y = d$PC2,
      z = d$PC3,
      type = "n",
      angle = 48,
      xlim = x_lim,
      ylim = y_lim,
      zlim = z_lim,
      xlab = "PC1 loading",
      ylab = "PC2 loading",
      zlab = "PC3 loading",
      main = species_labels[[sp]],
      grid = TRUE,
      box = FALSE
    )
    origin <- s3d$xyz.convert(0, 0, 0)
    points(origin$x, origin$y, pch = 16, cex = 0.8, col = "black")
    for (i in seq_len(nrow(d))) {
      end <- s3d$xyz.convert(d$PC1[i], d$PC2[i], d$PC3[i])
      arrows(
        origin$x, origin$y, end$x, end$y,
        col = species_cols[[sp]],
        lwd = 1.2,
        length = 0.07
      )
      text(end$x, end$y, labels = d$Variable[i], cex = 0.55, pos = 3, col = "grey20")
    }
  }

  plot.new()
  legend(
    "center",
    legend = c(
      "Arrow = BIO variable loading vector",
      paste0("Top ", top_n_arrows, " variables by |PC1,PC2,PC3| per species")
    ),
    col = c("grey20", "grey20"),
    lty = c(1, NA),
    lwd = c(1.2, NA),
    bty = "n",
    cex = 0.95,
    title = "Interpretation"
  )
  mtext("SuppFig2: Species-wise WorldClim PCA Loading Vectors (PC1, PC2, PC3)", outer = TRUE, cex = 1.2, font = 2)
  dev.off()
}

out2_png <- file.path(fig_dir, "SuppFig2_WorldClim_PCA_species_3D.png")
out2_pdf <- file.path(fig_dir, "SuppFig2_WorldClim_PCA_species_3D.pdf")
draw_suppfig2(out2_png, "png")
draw_suppfig2(out2_pdf, "pdf")

message("Saved: ", out1_png)
message("Saved: ", out1_pdf)
message("Saved: ", out2_png)
message("Saved: ", out2_pdf)
