suppressPackageStartupMessages({
  library(ggplot2)
})

base_dir <- '/Users/nirwantandukar/Documents/Github/Landadapt'
go_dir <- file.path(base_dir, 'results', 'tables', 'GO_terms')
fig_dir <- file.path(base_dir, 'Figs', 'Supplementary')
out_dir <- file.path(base_dir, 'results', 'tables', 'GO_terms')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

files <- data.frame(
  phenotype = c('PCs', 'pH', 'SoilN', 'AMF', 'Aridity'),
  BP = c('PCs_Maize_FDR_GO_BP.txt', 'pH_Maize_FDR_GO_BP.txt', 'SoilN_Maize_FDR_GO_BP.txt', 'AMF_Maize_FDR_GO_BP.txt', NA),
  MF = c('PCs_Maize_FDR_GO_MF.txt', 'pH_Maize_FDR_GO_MF.txt', 'SoilN_Maize_FDR_GO_MF.txt', 'AMF_Maize_FDR_GO_MF.txt', 'Aridity_Maize_FDR_GO_MF.txt'),
  stringsAsFactors = FALSE
)

pretty_pheno <- c(PCs = 'WorldClim PCs', pH = 'soil pH', SoilN = 'soil nitrogen', AMF = 'AM fungal phenotypes', Aridity = 'aridity index')
pretty_onto <- c(BP = 'GO biological process', MF = 'GO molecular function')
panel_label <- list(
  BP = c(PCs = 'A', pH = 'B', SoilN = 'C', AMF = 'D'),
  MF = c(PCs = 'A', pH = 'B', SoilN = 'C', AMF = 'D', Aridity = 'E')
)

wrap_term <- function(x, width = 38) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = '\n'), character(1))
}

reorder_within <- function(x, by, within, fun = mean, sep = '___', ...) {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)
}

scale_y_reordered <- function(..., sep = '___') {
  scale_y_discrete(labels = function(x) gsub(paste0(sep, '.*$'), '', x), ...)
}

parse_panther <- function(path, phenotype, ontology) {
  lines <- readLines(path, warn = FALSE)
  start_idx <- grep('^GO biological process complete|^GO molecular function complete', lines)
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
  if (ncol(tab) < 8) stop(sprintf('Unexpected column count in %s', path))
  tab <- tab[, 1:8]
  names(tab) <- c('term', 'ref_count', 'input_count', 'expected', 'direction', 'fold_enrichment', 'raw_p', 'fdr')
  tab$ref_count <- suppressWarnings(as.numeric(tab$ref_count))
  tab$input_count <- suppressWarnings(as.numeric(tab$input_count))
  tab$expected <- suppressWarnings(as.numeric(tab$expected))
  tab$fold_enrichment <- suppressWarnings(as.numeric(tab$fold_enrichment))
  tab$raw_p <- suppressWarnings(as.numeric(tab$raw_p))
  tab$fdr <- suppressWarnings(as.numeric(tab$fdr))
  tab <- subset(tab, !is.na(fdr) & direction == '+' & fdr < 0.05)
  if (nrow(tab) == 0) return(tab)
  tab$phenotype <- phenotype
  tab$ontology <- ontology
  tab
}

load_ontology <- function(ontology = c('BP', 'MF')) {
  ontology <- match.arg(ontology)
  out <- list()
  for (i in seq_len(nrow(files))) {
    f <- files[i, ontology]
    if (is.na(f) || f == '') next
    path <- file.path(go_dir, f)
    if (!file.exists(path)) next
    out[[length(out) + 1]] <- parse_panther(path, files$phenotype[i], ontology)
  }
  if (length(out) == 0) return(data.frame())
  do.call(rbind, out)
}

make_top_table <- function(df, n = 10) {
  pieces <- lapply(split(df, df$phenotype), function(d) {
    d <- d[order(-d$fold_enrichment, d$fdr, d$raw_p), ]
    head(d, n)
  })
  out <- do.call(rbind, pieces)
  out$phenotype_label <- pretty_pheno[out$phenotype]
  out$term_wrapped <- wrap_term(out$term, width = 34)
  out$term_facet <- reorder_within(out$term_wrapped, out$fold_enrichment, out$phenotype_label)
  out$neglog10_fdr <- -log10(out$fdr)
  out
}

make_plot <- function(df, ontology = c('BP', 'MF')) {
  ontology <- match.arg(ontology)
  if (nrow(df) == 0) return(NULL)
  df$facet_title <- paste0(panel_label[[ontology]][df$phenotype], '. ', pretty_pheno[df$phenotype])
  ncol_facets <- if (ontology == 'BP') 2 else 2
  ggplot(df, aes(x = fold_enrichment, y = term_facet)) +
    geom_segment(aes(x = 0, xend = fold_enrichment, yend = term_facet), color = 'grey75', linewidth = 0.4) +
    geom_point(aes(size = input_count, color = neglog10_fdr), alpha = 0.95) +
    facet_wrap(~facet_title, scales = 'free_y', ncol = ncol_facets) +
    scale_y_reordered() +
    scale_color_gradient(low = '#fdd49e', high = '#b30000') +
    scale_size(range = c(2.5, 8)) +
    labs(
      title = paste('Maize GO enrichment from shared ortholog genes in at least 3 species:', pretty_onto[[ontology]]),
      subtitle = 'Top 10 FDR-significant overrepresented terms per phenotype ranked by fold enrichment',
      x = 'Fold enrichment',
      y = NULL,
      color = '-log10(FDR)',
      size = 'Input\ngenes'
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = 'bold', size = 16),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(face = 'bold', size = 11),
      axis.text.y = element_text(size = 9),
      legend.position = 'right',
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

bp_df <- load_ontology('BP')
mf_df <- load_ontology('MF')

bp_top <- make_top_table(bp_df, n = 10)
mf_top <- make_top_table(mf_df, n = 10)

write.csv(bp_top, file.path(out_dir, 'Maize_shared_ge3_GO_BP_top10_fold_for_plot.csv'), row.names = FALSE)
write.csv(mf_top, file.path(out_dir, 'Maize_shared_ge3_GO_MF_top10_fold_for_plot.csv'), row.names = FALSE)

bp_plot <- make_plot(bp_top, 'BP')
mf_plot <- make_plot(mf_top, 'MF')

ggsave(file.path(fig_dir, 'SuppFig_GO_maize_shared_ge3_BP_top10_fold.png'), bp_plot, width = 14, height = 11, dpi = 300)
ggsave(file.path(fig_dir, 'SuppFig_GO_maize_shared_ge3_BP_top10_fold.pdf'), bp_plot, width = 14, height = 11)
ggsave(file.path(fig_dir, 'SuppFig_GO_maize_shared_ge3_MF_top10_fold.png'), mf_plot, width = 14, height = 13, dpi = 300)
ggsave(file.path(fig_dir, 'SuppFig_GO_maize_shared_ge3_MF_top10_fold.pdf'), mf_plot, width = 14, height = 13)

cat('WROTE', file.path(fig_dir, 'SuppFig_GO_maize_shared_ge3_BP_top10_fold.png'), '\n')
cat('WROTE', file.path(fig_dir, 'SuppFig_GO_maize_shared_ge3_MF_top10_fold.png'), '\n')
cat('WROTE', file.path(out_dir, 'Maize_shared_ge3_GO_BP_top10_fold_for_plot.csv'), '\n')
cat('WROTE', file.path(out_dir, 'Maize_shared_ge3_GO_MF_top10_fold_for_plot.csv'), '\n')
