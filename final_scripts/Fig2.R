suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

base_dir <- '/Users/nirwantandukar/Documents/Github/Landadapt'
summary_path <- file.path(base_dir, 'results', 'tables', 'Null_model', 'orthogroup_null_permutation_summary.csv')
dist_path <- file.path(base_dir, 'results', 'tables', 'Null_model', 'orthogroup_null_permutation_distributions.csv')
out_dir <- file.path(base_dir, 'Figs', 'Fig2')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- read.csv(summary_path, stringsAsFactors = FALSE)
dist_df <- read.csv(dist_path, stringsAsFactors = FALSE)

trait_order <- c('PC1', 'PC2', 'PC3', 'pH', 'soilN', 'AM_rel', 'AM_roots', 'aridity')
trait_labels <- c(PC1 = 'PC1', PC2 = 'PC2', PC3 = 'PC3', pH = 'pH', soilN = 'SoilN', AM_rel = 'AM rel', AM_roots = 'AM roots', aridity = 'Aridity')
palette_vals <- c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3')

summary_df$Trait_label <- factor(summary_df$Trait_label, levels = trait_order)
dist_df$Trait_label <- factor(dist_df$Trait_label, levels = trait_order)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(
      size = 14,
      face = 'bold',
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x = element_text(
      size = 24,
      face = 'bold'
    ),
    axis.title.y = element_text(
      size = 24,
      face = 'bold'
    ),
    axis.text.x = element_text(
      size = 24,
      color = 'black'
    ),
    axis.text.y = element_text(
      size = 24,
      color = 'black'
    ),
    axis.line = element_line(color = 'black'),
    panel.grid = element_blank(),
    legend.position = 'top',
    legend.title = element_blank(),
    legend.text = element_text(
      size = 16
    ),
    plot.margin = margin(15, 15, 15, 15)
  )

make_panel <- function(dist_col, obs_col, mean_col, low_col, high_col, p_col, title_text, show_legend = TRUE) {
  dat <- summary_df
  dat$obs_label <- sprintf('obs=%d\nP=%.4f', dat[[obs_col]], dat[[p_col]])
  dat$label_y <- pmax(dat[[obs_col]], dat[[high_col]]) + ifelse(dist_col == 'all5_count', 0.18, 0.6)

  p <- ggplot() +
    geom_violin(
      data = dist_df,
      aes(x = Trait_label, y = .data[[dist_col]], fill = Trait_label),
      color = 'black', linewidth = 0.4, trim = TRUE, alpha = 0.75,
      show.legend = FALSE
    ) +
    geom_errorbar(
      data = dat,
      aes(x = Trait_label, ymin = .data[[low_col]], ymax = .data[[high_col]], color = 'Null mean with 95% interval'),
      width = 0.08, linewidth = 0.7
    ) +
    geom_point(
      data = dat,
      aes(x = Trait_label, y = .data[[mean_col]], color = 'Null mean with 95% interval'),
      size = 2.8
    ) +
    geom_point(
      data = dat,
      aes(x = Trait_label, y = .data[[obs_col]], shape = 'Observed'),
      size = 4.5, color = 'black', fill = '#DC143C', stroke = 0.4
    ) +
    geom_text(
      data = dat,
      aes(x = Trait_label, y = label_y, label = obs_label),
      size = 4.4, vjust = 0
    ) +
    scale_fill_manual(values = setNames(palette_vals, trait_order), guide = 'none') +
    scale_shape_manual(values = c(Observed = 23), breaks = 'Observed') +
    scale_color_manual(values = c('Null mean with 95% interval' = 'black')) +
    scale_x_discrete(labels = trait_labels) +
    labs(title = title_text, x = NULL, y = 'Orthogroup count') +
    plot_theme +
    theme(
      plot.title = element_text(hjust = 0, face = 'bold', size = 18),
      axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.y = element_text(size = 18, face = 'bold'),
      legend.position = if (show_legend) 'top' else 'none',
      legend.background = element_rect(fill = alpha('white', 0.9), color = 'grey70'),
      legend.key = element_blank()
    ) +
    guides(
      shape = guide_legend(override.aes = list(fill = '#DC143C', color = 'black', size = 4)),
      color = guide_legend(override.aes = list(shape = 16, size = 3))
    )

  p
}

p1 <- make_panel('all5_count', 'Observed_all5', 'Null_mean_all5', 'Null_q025_all5', 'Null_q975_all5', 'Empirical_p_all5', 'A. All-five shared orthogroup recurrence vs null', TRUE) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p2 <- make_panel('ge4_count', 'Observed_ge4', 'Null_mean_ge4', 'Null_q025_ge4', 'Null_q975_ge4', 'Empirical_p_ge4', 'B. Shared orthogroups in >=4 species vs null', FALSE) +
  labs(x = 'Trait')

title_grob <- textGrob('Observed shared orthogroup recurrence exceeds matched-null expectations',
                       gp = gpar(fontsize = 20, fontface = 'bold'))
footer_grob <- textGrob('AM rel = AM fungal relative abundance colonization; AM roots = AM fungal roots colonized',
                        gp = gpar(fontsize = 12))

combined <- arrangeGrob(
  p1, p2,
  ncol = 1,
  heights = c(1, 1.2),
  top = title_grob,
  bottom = footer_grob
)

png_path <- file.path(out_dir, 'Fig2_null_permutation_orthogroup_recurrence.png')
pdf_path <- file.path(out_dir, 'Fig2_null_permutation_orthogroup_recurrence.pdf')

ggsave(png_path, combined, width = 14, height = 12, dpi = 300, bg = 'white')
ggsave(pdf_path, combined, width = 14, height = 12, bg = 'white')

cat('WROTE', png_path, '\n')
cat('WROTE', pdf_path, '\n')
