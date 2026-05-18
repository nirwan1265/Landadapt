suppressPackageStartupMessages({
  library(ggplot2)
})

base_dir <- '/Users/nirwantandukar/Documents/Github/Landadapt'
go_dir <- file.path(base_dir, 'results', 'tables', 'GO_terms')
out_dir <- file.path(base_dir, 'Figs', 'Fig4')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- data.frame(
  phenotype = c('PCs', 'pH', 'SoilN', 'AMF'),
  BP = c('PCs_Maize_FDR_GO_BP.txt', 'pH_Maize_FDR_GO_BP.txt', 'SoilN_Maize_FDR_GO_BP.txt', 'AMF_Maize_FDR_GO_BP.txt'),
  stringsAsFactors = FALSE
)

pretty_pheno <- c(PCs = 'WorldClim PCs', pH = 'soil pH', SoilN = 'soil nitrogen', AMF = 'AM fungal phenotypes')
panel_label <- c(PCs = 'A', pH = 'B', SoilN = 'C', AMF = 'D')

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

wrap_term <- function(x, width = 34) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = '\n'), character(1))
}

reorder_within <- function(x, by, within, fun = mean, sep = '___', ...) {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)
}

scale_y_reordered <- function(..., sep = '___') {
  scale_y_discrete(labels = function(x) gsub(paste0(sep, '.*$'), '', x), ...)
}

parse_panther <- function(path, phenotype) {
  lines <- readLines(path, warn = FALSE)
  start_idx <- grep('^GO biological process complete', lines)
  if (length(start_idx) == 0) stop(sprintf('Could not find header row in %s', path))
  tab <- read.delim(
    text = paste(lines[(start_idx[1] + 1):length(lines)], collapse = '\n'),
    header = FALSE,
    sep = '\t',
    quote = '',
    fill = TRUE,
    stringsAsFactors = FALSE,
    comment.char = ''
  )
  tab <- tab[, 1:8]
  names(tab) <- c('term', 'ref_count', 'input_count', 'expected', 'direction', 'fold_enrichment', 'raw_p', 'fdr')
  tab$input_count <- suppressWarnings(as.numeric(tab$input_count))
  tab$fold_enrichment <- suppressWarnings(as.numeric(tab$fold_enrichment))
  tab$raw_p <- suppressWarnings(as.numeric(tab$raw_p))
  tab$fdr <- suppressWarnings(as.numeric(tab$fdr))
  tab <- subset(tab, !is.na(fdr) & direction == '+' & fdr < 0.05)
  tab$phenotype <- phenotype
  tab
}

all_bp <- do.call(rbind, lapply(seq_len(nrow(files)), function(i) {
  parse_panther(file.path(go_dir, files$BP[i]), files$phenotype[i])
}))

bp_top <- do.call(rbind, lapply(split(all_bp, all_bp$phenotype), function(d) {
  d <- d[order(-d$fold_enrichment, d$fdr, d$raw_p), ]
  head(d, 10)
}))

bp_top$phenotype_label <- pretty_pheno[bp_top$phenotype]
bp_top$facet_title <- paste0(panel_label[bp_top$phenotype], '. ', pretty_pheno[bp_top$phenotype])
bp_top$term_wrapped <- wrap_term(bp_top$term, width = 34)
bp_top$term_facet <- reorder_within(bp_top$term_wrapped, bp_top$fold_enrichment, bp_top$phenotype_label)
bp_top$neglog10_fdr <- -log10(bp_top$fdr)

write.csv(bp_top, file.path(go_dir, 'Maize_shared_ge3_GO_BP_top10_fold_for_plot.csv'), row.names = FALSE)

p <- ggplot(bp_top, aes(x = fold_enrichment, y = term_facet)) +
  geom_segment(aes(x = 0, xend = fold_enrichment, yend = term_facet), color = 'grey75', linewidth = 0.4) +
  geom_point(aes(size = input_count, color = neglog10_fdr), alpha = 0.95) +
  facet_wrap(~facet_title, scales = 'free_y', ncol = 2) +
  scale_y_reordered() +
  scale_color_gradient(low = '#fdd49e', high = '#b30000') +
  scale_size(range = c(2.5, 8)) +
  labs(
    title = 'Maize GO enrichment from shared ortholog genes in at least 3 species: GO biological process',
    subtitle = 'Top 10 FDR-significant overrepresented terms per phenotype ranked by fold enrichment',
    x = 'Fold enrichment',
    y = NULL,
    color = '-log10(FDR)',
    size = 'Input\ngenes'
  ) +
  plot_theme +
  theme(
    plot.title = element_text(face = 'bold', size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    strip.text = element_text(face = 'bold', size = 11),
    axis.text.y = element_text(size = 9, color = 'black'),
    axis.text.x = element_text(size = 12, color = 'black'),
    axis.title.x = element_text(size = 14, face = 'bold'),
    legend.position = 'right',
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

png_path <- file.path(out_dir, 'Fig4_GO_maize_shared_ge3_BP_top10_fold.png')
pdf_path <- file.path(out_dir, 'Fig4_GO_maize_shared_ge3_BP_top10_fold.pdf')

ggsave(png_path, p, width = 14, height = 11, dpi = 300, bg = 'white')
ggsave(pdf_path, p, width = 14, height = 11, bg = 'white')

cat('WROTE', png_path, '\n')
cat('WROTE', pdf_path, '\n')
