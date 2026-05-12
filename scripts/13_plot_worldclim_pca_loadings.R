suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(stringr)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
loadings_dir <- "data/WorldClim"
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

write_csv(bio_names, file.path(res_dir, "worldclim_bio_variable_names.csv"))

# ------------------------------------------------------------
# Load and reshape species loadings
# ------------------------------------------------------------
files <- list.files(loadings_dir, pattern = "^worldclim_.*_PCA_loadings\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No PCA loading files found in ", loadings_dir)

load_one <- function(path) {
  species <- basename(path) %>%
    str_remove("^worldclim_") %>%
    str_remove("_PCA_loadings\\.csv$")

  read_csv(path, show_col_types = FALSE) %>%
    pivot_longer(cols = c(PC1, PC2, PC3), names_to = "PC", values_to = "loading") %>%
    mutate(species = species)
}

loadings_long <- bind_rows(lapply(files, load_one)) %>%
  left_join(bio_names, by = "Variable") %>%
  mutate(
    species = factor(species, levels = c("arabidopsis", "barley", "rice", "maize", "sorghum")),
    PC = factor(PC, levels = c("PC1", "PC2", "PC3")),
    bio_label = paste0(Variable, "  ", bio_name)
  )

write_csv(loadings_long, file.path(res_dir, "worldclim_pca_loadings_long.csv"))

# ------------------------------------------------------------
# Figure 1: full loading heatmap by species
# ------------------------------------------------------------
lim <- max(abs(loadings_long$loading), na.rm = TRUE)

p_species <- ggplot(loadings_long, aes(x = PC, y = fct_rev(Variable), fill = loading)) +
  geom_tile(color = "grey95", linewidth = 0.25) +
  facet_wrap(~ species, ncol = 3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-lim, lim), name = "Loading"
  ) +
  labs(
    title = "WorldClim PCA Loadings by Species",
    x = NULL,
    y = "BIO variable"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(fig_dir, "WorldClim_PCA_loadings_by_species_heatmap.png"),
  p_species, width = 12, height = 9, dpi = 300, bg = "white", device = ragg::agg_png
)

# ------------------------------------------------------------
# Aggregate summaries across species
# ------------------------------------------------------------
summary_mean_abs <- loadings_long %>%
  group_by(PC, Variable, bio_name) %>%
  summarise(
    mean_abs_loading = mean(abs(loading), na.rm = TRUE),
    median_abs_loading = median(abs(loading), na.rm = TRUE),
    mean_signed_loading = mean(loading, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_mean_abs, file.path(res_dir, "worldclim_pc_loading_summary_meanAbs.csv"))

order_tbl <- summary_mean_abs %>%
  group_by(Variable, bio_name) %>%
  summarise(overall_mean_abs = mean(mean_abs_loading), .groups = "drop") %>%
  arrange(desc(overall_mean_abs))

summary_plot_df <- summary_mean_abs %>%
  left_join(order_tbl, by = c("Variable", "bio_name")) %>%
  mutate(bio_label = factor(paste0(Variable, "  ", bio_name), levels = rev(paste0(order_tbl$Variable, "  ", order_tbl$bio_name))))

# PC-specific ranking labels so each PC panel is independently ordered high -> low.
summary_pc_ranked <- summary_mean_abs %>%
  mutate(bio_label = paste0(Variable, "  ", bio_name)) %>%
  group_by(PC) %>%
  arrange(desc(mean_abs_loading), .by_group = TRUE) %>%
  mutate(bio_label_pc = factor(paste0(PC, "__", bio_label), levels = rev(unique(paste0(PC, "__", bio_label))))) %>%
  ungroup()

# ------------------------------------------------------------
# Figure 2: species-wise |loading| heatmap (no cross-species mean)
# ------------------------------------------------------------
species_abs_df <- loadings_long %>%
  mutate(
    abs_loading = abs(loading),
    bio_label = paste0(Variable, "  ", bio_name)
  ) %>%
  group_by(PC, Variable, bio_name, bio_label) %>%
  mutate(pc_rank = mean(abs_loading, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(PC) %>%
  arrange(desc(pc_rank), .by_group = TRUE) %>%
  mutate(bio_label = factor(bio_label, levels = unique(bio_label))) %>%
  ungroup()

p_species_abs <- ggplot(species_abs_df, aes(x = bio_label, y = species, fill = abs_loading)) +
  geom_tile(color = "grey95", linewidth = 0.25) +
  facet_wrap(~ PC, ncol = 1) +
  scale_fill_gradient(low = "white", high = "black", name = "|Loading|") +
  labs(
    title = "WorldClim PCA Loadings by Species (Absolute Value)",
    x = "BIO variable",
    y = "Species"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(fig_dir, "WorldClim_PCA_loadings_species_rows_heatmap.png"),
  p_species_abs, width = 14, height = 9, dpi = 300, bg = "white", device = ragg::agg_png
)

# ------------------------------------------------------------
# Figure 2B: mean |loading| by PC in 3 side-by-side columns
# ------------------------------------------------------------
p_mean_abs <- ggplot(summary_pc_ranked, aes(x = 1, y = bio_label_pc, fill = mean_abs_loading)) +
  geom_tile(color = "grey95", linewidth = 0.25) +
  facet_wrap(~ PC, ncol = 3, scales = "free_y") +
  scale_fill_gradient(low = "white", high = "black", name = "Mean |loading|") +
  scale_y_discrete(labels = function(x) sub("^[A-Z0-9]+__", "", x)) +
  labs(
    title = "WorldClim PCA Mean Absolute Loadings by PC",
    x = NULL,
    y = "BIO variable"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(fig_dir, "WorldClim_PCA_meanAbs_loadings_heatmap.png"),
  p_mean_abs, width = 14, height = 9, dpi = 300, bg = "white", device = ragg::agg_png
)

# ------------------------------------------------------------
# Figure 3: top contributors per PC
# ------------------------------------------------------------
top_n <- 8
top_per_pc <- summary_mean_abs %>%
  group_by(PC) %>%
  arrange(desc(mean_abs_loading), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  mutate(
    label = paste0(Variable, "  ", bio_name),
    label = fct_reorder(label, mean_abs_loading)
  )

p_top <- ggplot(top_per_pc, aes(x = mean_abs_loading, y = label, fill = PC)) +
  geom_col(width = 0.7) +
  facet_wrap(~ PC, scales = "free_y") +
  labs(
    title = "Top WorldClim Contributors to Each PC",
    x = "Mean |loading| across species",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  file.path(fig_dir, "WorldClim_PCA_top_contributors_barplot.png"),
  p_top, width = 13, height = 6, dpi = 300, bg = "white", device = ragg::agg_png
)

message("Saved figures in: ", fig_dir)
message("Saved tables in: ", res_dir)
