suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

base_dir <- '/Users/nirwantandukar/Documents/Github/Landadapt'
values_path <- file.path(base_dir, 'results', 'tables', 'Supplementary', 'FigX_cross_species_summary_values.csv')
out_dir <- file.path(base_dir, 'Figs', 'Fig3')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vals <- read.csv(values_path, stringsAsFactors = FALSE)

species_order <- c('Arabidopsis', 'Barley', 'Rice', 'Maize', 'Sorghum')
trait_order <- c('PC1', 'PC2', 'PC3', 'pH', 'soilN', 'AM_rel', 'AM_roots', 'aridity')
trait_labels <- c(PC1 = 'PC1', PC2 = 'PC2', PC3 = 'PC3', pH = 'pH', soilN = 'SoilN', AM_rel = 'AM rel', AM_roots = 'AM roots', aridity = 'Aridity')
bar_colors <- c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3')

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

make_heat_df <- function(panel_code) {
  d <- subset(vals, panel == panel_code)
  d$species <- factor(d$species, levels = rev(species_order))
  d$trait <- factor(d$trait, levels = trait_order)
  d$plot_value <- log10(d$value + 1)
  d$label <- format(d$value, big.mark = ',', trim = TRUE, scientific = FALSE)
  d$label_color <- ifelse(d$plot_value >= (min(d$plot_value) + 0.55 * (max(d$plot_value) - min(d$plot_value))), 'white', 'black')
  d
}

make_bar_df <- function(panel_code) {
  d <- subset(vals, panel == panel_code)
  d$trait <- factor(d$trait, levels = trait_order)
  d$trait_label <- trait_labels[d$trait]
  d$fill_col <- bar_colors[match(as.character(d$trait), trait_order)]
  d
}

A <- make_heat_df('A')
B <- make_heat_df('B')
C <- make_bar_df('C')
D <- make_bar_df('D')

pA <- ggplot(A, aes(x = trait, y = species, fill = plot_value)) +
  geom_tile(color = 'white', linewidth = 0.5) +
  geom_text(aes(label = label, color = label_color), size = 4) +
  scale_color_identity() +
  scale_fill_gradient(low = '#ffffcc', high = '#b10026') +
  scale_x_discrete(labels = trait_labels) +
  labs(title = 'A. Bonferroni-significant SNP counts', x = 'Trait', y = 'Species', fill = 'log10(count + 1)') +
  plot_theme +
  theme(
    plot.title = element_text(hjust = 0, face = 'bold', size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    legend.position = 'right',
    panel.grid = element_blank()
  )

pB <- ggplot(B, aes(x = trait, y = species, fill = plot_value)) +
  geom_tile(color = 'white', linewidth = 0.5) +
  geom_text(aes(label = label, color = label_color), size = 4) +
  scale_color_identity() +
  scale_fill_gradient(low = '#ffffd9', high = '#081d58') +
  scale_x_discrete(labels = trait_labels) +
  labs(title = 'B. Top-0.5% mapped gene counts', x = 'Trait', y = 'Species', fill = 'log10(count + 1)') +
  plot_theme +
  theme(
    plot.title = element_text(hjust = 0, face = 'bold', size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    legend.position = 'right',
    panel.grid = element_blank()
  )

pC <- ggplot(C, aes(x = trait, y = value, fill = trait)) +
  geom_col(color = 'black', linewidth = 0.4, show.legend = FALSE) +
  geom_text(aes(label = value), vjust = -0.25, size = 4) +
  scale_fill_manual(values = setNames(bar_colors, trait_order)) +
  scale_x_discrete(labels = trait_labels) +
  labs(title = 'C. All-five shared orthogroups', x = NULL, y = 'Orthogroup count') +
  plot_theme +
  theme(
    plot.title = element_text(hjust = 0, face = 'bold', size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16, face = 'bold')
  )

pD <- ggplot(D, aes(x = trait, y = value, fill = trait)) +
  geom_col(color = 'black', linewidth = 0.4, show.legend = FALSE) +
  geom_text(aes(label = value), vjust = -0.25, size = 4) +
  scale_fill_manual(values = setNames(bar_colors, trait_order)) +
  scale_x_discrete(labels = trait_labels) +
  labs(title = 'D. Shared orthogroups in >=4 species', x = NULL, y = 'Orthogroup count') +
  plot_theme +
  theme(
    plot.title = element_text(hjust = 0, face = 'bold', size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16, face = 'bold')
  )

main_title <- textGrob('Cross-species environmental GWAS summary across traits', gp = gpar(fontsize = 20, fontface = 'bold'))
footer <- textGrob('AM rel = AM fungal relative abundance colonization; AM roots = AM fungal roots colonized', gp = gpar(fontsize = 12))

combined <- arrangeGrob(
  grobs = list(pA, pC, pB, pD),
  layout_matrix = rbind(c(1, 2), c(3, 4)),
  widths = c(1.35, 1),
  top = main_title,
  bottom = footer
)

png_path <- file.path(out_dir, 'Fig3_cross_species_gwas_summary.png')
pdf_path <- file.path(out_dir, 'Fig3_cross_species_gwas_summary.pdf')

ggsave(png_path, combined, width = 15, height = 11, dpi = 300, bg = 'white')
ggsave(pdf_path, combined, width = 15, height = 11, bg = 'white')

cat('WROTE', png_path, '\n')
cat('WROTE', pdf_path, '\n')
