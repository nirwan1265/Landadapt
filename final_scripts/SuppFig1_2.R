suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

base_dir <- '/Users/nirwantandukar/Documents/Github/Landadapt'
wc_dir <- file.path(base_dir, 'data', 'WorldClim')
out_dir <- file.path(base_dir, 'Figs', 'Supplementary')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_order <- c('arabidopsis', 'barley', 'rice', 'maize', 'sorghum')
species_pretty <- c(arabidopsis = 'Arabidopsis', barley = 'Barley', rice = 'Rice', maize = 'Maize', sorghum = 'Sorghum')
species_cols <- c(arabidopsis = '#1b9e77', barley = '#7570b3', rice = '#e7298a', maize = '#d95f02', sorghum = '#66a61e')

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(
      size = 14,
      face = 'bold',
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x = element_text(size = 24, face = 'bold'),
    axis.title.y = element_text(size = 24, face = 'bold'),
    axis.text.x = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 24, color = 'black'),
    axis.line = element_line(color = 'black'),
    panel.grid = element_blank(),
    legend.position = 'top',
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(15, 15, 15, 15)
  )

# ---------- SuppFig1 ----------
var_list <- lapply(species_order, function(sp) {
  path <- file.path(wc_dir, paste0('worldclim_', sp, '_PCA_variance.csv'))
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$species <- sp
  df$PC_num <- seq_len(nrow(df))
  df
})
var_df <- do.call(rbind, var_list)
var_df$VariancePct <- 100 * var_df$VarianceExplained

summary_var <- aggregate(VariancePct ~ PC_num, data = var_df, FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
summary_var <- do.call(data.frame, summary_var)
colnames(summary_var) <- c('PC_num', 'Mean', 'SE')
summary_var <- subset(summary_var, PC_num <= 10)
summary_var$PC <- paste0('PC', summary_var$PC_num)
summary_var$PC <- factor(summary_var$PC, levels = paste0('PC', 1:10))

load_list <- lapply(species_order, function(sp) {
  path <- file.path(wc_dir, paste0('worldclim_', sp, '_PCA_loadings.csv'))
  df <- read.csv(path, stringsAsFactors = FALSE)
  keep <- df[, c('Variable', 'PC1', 'PC2', 'PC3')]
  long <- rbind(
    data.frame(species = sp, Variable = keep$Variable, PC = 'PC1', Loading = abs(keep$PC1)),
    data.frame(species = sp, Variable = keep$Variable, PC = 'PC2', Loading = abs(keep$PC2)),
    data.frame(species = sp, Variable = keep$Variable, PC = 'PC3', Loading = abs(keep$PC3))
  )
  long
})
load_df <- do.call(rbind, load_list)
load_df$species <- factor(load_df$species, levels = rev(species_order), labels = rev(species_pretty))
load_df$PC <- factor(load_df$PC, levels = c('PC1', 'PC2', 'PC3'))
load_df$Variable <- factor(load_df$Variable, levels = paste0('BIO', sprintf('%02d', 1:19)))

p1a <- ggplot(summary_var, aes(x = PC, y = Mean)) +
  geom_col(fill = '#636363', width = 0.72) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.18, linewidth = 0.5) +
  geom_line(aes(group = 1), color = '#636363', linewidth = 0.6) +
  geom_point(size = 1.5, color = '#636363') +
  labs(title = 'A', x = NULL, y = 'Mean variance explained (%)') +
  plot_theme +
  theme(
    plot.title = element_text(hjust = 0, size = 12, face = 'plain'),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 10, face = 'bold'),
    axis.line = element_blank(),
    panel.grid.major.y = element_line(color = 'grey85', linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

p1b <- ggplot(load_df, aes(x = Variable, y = species, fill = Loading)) +
  geom_tile(color = 'grey80', linewidth = 0.25) +
  facet_wrap(~PC, ncol = 1) +
  scale_fill_gradient(low = 'white', high = 'grey20') +
  labs(title = 'B', x = 'BIO variable', y = 'Species', fill = '|Loading|') +
  plot_theme +
  theme(
    plot.title = element_text(hjust = 0, size = 12, face = 'plain'),
    axis.text.x = element_text(size = 6, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 9, face = 'bold'),
    axis.title.y = element_text(size = 9, face = 'bold'),
    strip.text = element_text(size = 9, face = 'bold'),
    legend.position = 'right'
  )

supp1 <- arrangeGrob(
  p1a, p1b,
  ncol = 1,
  heights = c(1, 1.25),
  top = textGrob('Shared PCA Scree Across Species and Absolute WorldClim Loadings', gp = gpar(fontsize = 16, fontface = 'bold'))
)

# ---------- SuppFig2 ----------
project_iso <- function(x, y, z) {
  list(
    x = x - 0.55 * y,
    y = z + 0.35 * y
  )
}

draw_species_panel <- function(sp) {
  path <- file.path(wc_dir, paste0('worldclim_', sp, '_PCA_loadings.csv'))
  df <- read.csv(path, stringsAsFactors = FALSE)
  ord <- order(rowSums(abs(df[, c('PC1', 'PC2', 'PC3')])), decreasing = TRUE)
  df <- df[ord[1:min(12, nrow(df))], ]
  pr <- project_iso(df$PC1, df$PC2, df$PC3)
  pts <- data.frame(label = df$Variable, x = pr$x, y = pr$y)
  origin <- project_iso(0,0,0)

  axis_vals <- seq(-0.8, 0.8, by = 0.4)
  pc1_axis <- data.frame(x = c(-0.8, 0.8), y = c(0, 0), z = c(0, 0))
  pc2_axis <- data.frame(x = c(0, 0), y = c(-0.8, 0.8), z = c(0, 0))
  pc3_axis <- data.frame(x = c(0, 0), y = c(0, 0), z = c(-0.8, 0.8))
  pc1_proj <- project_iso(pc1_axis$x, pc1_axis$y, pc1_axis$z)
  pc2_proj <- project_iso(pc2_axis$x, pc2_axis$y, pc2_axis$z)
  pc3_proj <- project_iso(pc3_axis$x, pc3_axis$y, pc3_axis$z)

  pc1_ticks <- data.frame(val = axis_vals, lab = sprintf('%.1f', axis_vals))
  pc1_tick_pos <- project_iso(pc1_ticks$val, 0, 0)
  pc1_ticks$x <- pc1_tick_pos$x
  pc1_ticks$y <- pc1_tick_pos$y

  pc2_ticks <- data.frame(val = axis_vals, lab = sprintf('%.1f', axis_vals))
  pc2_tick_pos <- project_iso(0, pc2_ticks$val, 0)
  pc2_ticks$x <- pc2_tick_pos$x
  pc2_ticks$y <- pc2_tick_pos$y

  pc3_ticks <- data.frame(val = axis_vals, lab = sprintf('%.1f', axis_vals))
  pc3_tick_pos <- project_iso(0, 0, pc3_ticks$val)
  pc3_ticks$x <- pc3_tick_pos$x
  pc3_ticks$y <- pc3_tick_pos$y

  ggplot() +
    annotate('segment', x = pc1_proj$x[1], y = pc1_proj$y[1], xend = pc1_proj$x[2], yend = pc1_proj$y[2], color = 'grey35', linewidth = 0.45) +
    annotate('segment', x = pc2_proj$x[1], y = pc2_proj$y[1], xend = pc2_proj$x[2], yend = pc2_proj$y[2], color = 'grey35', linewidth = 0.45) +
    annotate('segment', x = pc3_proj$x[1], y = pc3_proj$y[1], xend = pc3_proj$x[2], yend = pc3_proj$y[2], color = 'grey35', linewidth = 0.45) +
    geom_segment(data = pts, aes(x = origin$x, y = origin$y, xend = x, yend = y),
                 arrow = arrow(length = unit(0.12, 'cm')), color = species_cols[[sp]], linewidth = 0.5) +
    geom_text(data = pts, aes(x = x, y = y, label = label), size = 2.2, color = 'grey30', vjust = -0.3) +
    geom_text(data = pc1_ticks, aes(x = x, y = y - 0.06, label = lab), size = 2.1, color = 'grey30') +
    geom_text(data = pc2_ticks, aes(x = x + 0.07, y = y, label = lab), size = 2.1, color = 'grey30') +
    geom_text(data = pc3_ticks, aes(x = x - 0.08, y = y, label = lab), size = 2.1, color = 'grey30') +
    annotate('text', x = 0, y = -1.02, label = 'PC1 loading', size = 3) +
    annotate('text', x = -1.02, y = 0.05, label = 'PC3 loading', angle = 90, size = 3) +
    annotate('text', x = 1.02, y = 0.08, label = 'PC2 loading', angle = 90, size = 3) +
    ggtitle(species_pretty[[sp]]) +
    coord_equal(xlim = c(-1.12, 1.12), ylim = c(-1.08, 0.95), expand = FALSE) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = 'bold'),
      plot.margin = margin(16, 18, 16, 18)
    )
}

p_list <- lapply(species_order, draw_species_panel)
legend_grob <- textGrob('Interpretation\nArrow = variable loading vector\nTop 12 variables by |PC1|, |PC2|, |PC3| per species',
                        x = 0.05, y = 0.8, hjust = 0, gp = gpar(fontsize = 10))
blank <- rectGrob(gp = gpar(col = NA, fill = 'white'))
layout <- rbind(c(1,2,3), c(4,5,6))
# convert legend_grob into grob slot 6 by making arranged list
supp2 <- arrangeGrob(
  grobs = c(p_list, list(arrangeGrob(blank, top = legend_grob))),
  layout_matrix = layout,
  top = textGrob('SuppFig2: Species-wise WorldClim PCA Loading Vectors (PC1, PC2, PC3)', gp = gpar(fontsize = 16, fontface = 'bold'))
)

png1 <- file.path(out_dir, 'SuppFig1_WorldClim_PCA.png')
pdf1 <- file.path(out_dir, 'SuppFig1_WorldClim_PCA.pdf')
png2 <- file.path(out_dir, 'SuppFig2_WorldClim_PCA_species_3D.png')
pdf2 <- file.path(out_dir, 'SuppFig2_WorldClim_PCA_species_3D.pdf')

ggsave(png1, supp1, width = 8, height = 8, dpi = 300, bg = 'white')
ggsave(pdf1, supp1, width = 8, height = 8, bg = 'white')
ggsave(png2, supp2, width = 12, height = 8, dpi = 300, bg = 'white')
ggsave(pdf2, supp2, width = 12, height = 8, bg = 'white')

cat('WROTE', png1, '\n')
cat('WROTE', png2, '\n')
