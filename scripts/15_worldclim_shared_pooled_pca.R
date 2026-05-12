suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(patchwork)
  library(stringr)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
in_dir <- "data/WorldClim"
fig_dir <- "Figs/Supplementary"
res_dir <- "results"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# BIO variable names
# ------------------------------------------------------------
bio_names <- tibble(
  Variable = sprintf("BIO%02d", 1:19),
  bio_name = c(
    "Annual Mean Temperature",
    "Mean Diurnal Range",
    "Isothermality",
    "Temperature Seasonality",
    "Max Temperature of Warmest Month",
    "Min Temperature of Coldest Month",
    "Temperature Annual Range",
    "Mean Temperature of Wettest Quarter",
    "Mean Temperature of Driest Quarter",
    "Mean Temperature of Warmest Quarter",
    "Mean Temperature of Coldest Quarter",
    "Annual Precipitation",
    "Precipitation of Wettest Month",
    "Precipitation of Driest Month",
    "Precipitation Seasonality",
    "Precipitation of Wettest Quarter",
    "Precipitation of Driest Quarter",
    "Precipitation of Warmest Quarter",
    "Precipitation of Coldest Quarter"
  )
)

# ------------------------------------------------------------
# Load species files and pool
# ------------------------------------------------------------
files <- list.files(in_dir, pattern = "^worldclim_.*_BIO_PCA\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No worldclim_*_BIO_PCA.csv files found in ", in_dir)

read_one <- function(path) {
  sp <- basename(path) %>%
    str_remove("^worldclim_") %>%
    str_remove("_BIO_PCA\\.csv$")

  dat <- read_csv(path, show_col_types = FALSE)
  req <- sprintf("BIO%02d", 1:19)
  if (!all(req %in% names(dat))) {
    stop("Missing BIO01..BIO19 columns in: ", path)
  }

  dat %>%
    mutate(
      species = sp,
      Lines = as.character(Lines),
      Long = as.numeric(Long),
      Lat = as.numeric(Lat)
    ) %>%
    select(species, Lines, Long, Lat, all_of(req))
}

pooled <- bind_rows(lapply(files, read_one))
bio_cols <- sprintf("BIO%02d", 1:19)
pooled_clean <- pooled %>%
  filter(if_all(all_of(bio_cols), ~ !is.na(.x)))

if (nrow(pooled_clean) < 10) {
  stop("Too few complete rows after NA filtering for pooled PCA.")
}

# ------------------------------------------------------------
# Shared pooled PCA
# ------------------------------------------------------------
bio_mat <- as.matrix(pooled_clean[, bio_cols])
pca <- prcomp(bio_mat, center = TRUE, scale. = TRUE)

var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as_tibble(pca$x) %>%
  bind_cols(pooled_clean %>% select(species, Lines, Long, Lat))

loadings <- as_tibble(pca$rotation, rownames = "Variable")

write_csv(
  scores,
  file.path(res_dir, "worldclim_shared_pooled_pca_scores.csv")
)
write_csv(
  loadings,
  file.path(res_dir, "worldclim_shared_pooled_pca_loadings.csv")
)
write_csv(
  tibble(PC = paste0("PC", seq_along(var_expl)), var_explained = var_expl),
  file.path(res_dir, "worldclim_shared_pooled_pca_variance.csv")
)

# ------------------------------------------------------------
# Figure data
# ------------------------------------------------------------
load_long <- loadings %>%
  select(Variable, PC1, PC2, PC3) %>%
  pivot_longer(cols = c(PC1, PC2, PC3), names_to = "PC", values_to = "loading") %>%
  left_join(bio_names, by = "Variable") %>%
  mutate(
    label = paste0(Variable, "  ", bio_name),
    label = factor(label, levels = rev(paste0(sprintf("BIO%02d", 1:19), "  ", bio_names$bio_name))),
    PC = factor(
      PC,
      levels = c("PC1", "PC2", "PC3"),
      labels = c(
        paste0("PC1 (", sprintf("%.1f", 100 * var_expl[1]), "%)"),
        paste0("PC2 (", sprintf("%.1f", 100 * var_expl[2]), "%)"),
        paste0("PC3 (", sprintf("%.1f", 100 * var_expl[3]), "%)")
      )
    )
  )

scree_df <- tibble(
  PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  var_explained = 100 * var_expl[1:10]
)

species_levels <- c("arabidopsis", "barley", "maize", "rice", "sorghum")
scores_plot <- scores %>%
  mutate(species = factor(species, levels = species_levels))

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------
p_scree <- ggplot(scree_df, aes(x = PC, y = var_explained, group = 1)) +
  geom_col(fill = "grey35", width = 0.75) +
  geom_line(color = "black") +
  geom_point(color = "black", size = 1.8) +
  labs(title = "A. Shared PCA Scree", x = NULL, y = "Variance explained (%)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

lim <- max(abs(load_long$loading), na.rm = TRUE)
p_load <- ggplot(load_long, aes(x = PC, y = label, fill = loading)) +
  geom_tile(color = "grey95", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-lim, lim), name = "Loading"
  ) +
  labs(title = "B. BIO Loadings in Shared PCA", x = NULL, y = "BIO variable") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

p_species <- ggplot(scores_plot, aes(x = species, y = PC1, fill = species)) +
  geom_violin(scale = "width", color = "grey20", linewidth = 0.25, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.8) +
  scale_fill_manual(values = c(
    arabidopsis = "#009E73",
    barley = "#0072B2",
    maize = "#E69F00",
    rice = "#CC79A7",
    sorghum = "#D55E00"
  )) +
  labs(
    title = "C. Shared PC1 Score Distribution by Species",
    x = NULL,
    y = "Shared PC1 score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = "bold")
  )

fig <- (p_scree | p_load) / p_species +
  plot_annotation(
    title = "Shared Pooled WorldClim PCA Across All Species",
    subtitle = "PCA fit on pooled BIO01-BIO19 values; common axis definitions across species"
  )

out_fig <- file.path(fig_dir, "Supp_shared_pooled_WorldClim_PCA.png")
ggsave(out_fig, fig, width = 14, height = 11, dpi = 300, bg = "white", device = ragg::agg_png)

message("Saved: ", out_fig)
message("Saved: ", file.path(res_dir, "worldclim_shared_pooled_pca_scores.csv"))
message("Saved: ", file.path(res_dir, "worldclim_shared_pooled_pca_loadings.csv"))
message("Saved: ", file.path(res_dir, "worldclim_shared_pooled_pca_variance.csv"))
