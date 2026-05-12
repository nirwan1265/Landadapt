suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(data.table)
})

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------
fig_dir <- "Figs/Fig4"
res_dir <- "results/pc1_worldclim"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

top_fraction <- 0.005
window_bp <- 25000L

species_order <- c("Arabidopsis", "Barley", "Rice", "Maize", "Sorghum")
species_to_crop <- c(
  Arabidopsis = "arabidopsis",
  Barley = "barley",
  Rice = "rice",
  Maize = "maize",
  Sorghum = "sorghum"
)

crop_colors <- c(
  maize = "#E69F00",
  sorghum = "#D55E00",
  rice = "#CC79A7",
  arabidopsis = "#009E73",
  barley = "#0072B2",
  pearl_millet = "#F0E442"
)

gwas_cfg <- list(
  list(species = "Arabidopsis", file = "/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/arabidopsis/arabidopsis_PC1_worldclim.txt", sep = "\t", chr_col = "chr", pos_col = "ps", p_col = "p_wald", snp_col = "rs"),
  list(species = "Barley", file = "/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/barley/barley_PC1_worldclim.txt", sep = "\t", chr_col = "chr", pos_col = "ps", p_col = "p_wald", snp_col = "rs"),
  list(species = "Rice", file = "/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/rice/rice_PC1_worldclim.txt", sep = "\t", chr_col = "chr", pos_col = "ps", p_col = "p_wald", snp_col = "rs"),
  list(species = "Maize", file = "/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/maize/maize_PC1_worldclim.txt", sep = "\t", chr_col = "chr", pos_col = "ps", p_col = "p_wald", snp_col = "rs"),
  list(species = "Sorghum", file = "/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/sorghum/sorghum_PC1_worldclim.txt", sep = "\t", chr_col = "chr", pos_col = "ps", p_col = "p_wald", snp_col = "rs")
)

pheno_cfg <- list(
  list(species = "Arabidopsis", file = "data/WorldClim/worldclim_arabidopsis_BIO_PCA.csv"),
  list(species = "Barley", file = "data/WorldClim/worldclim_barley_BIO_PCA.csv"),
  list(species = "Rice", file = "data/WorldClim/worldclim_rice_BIO_PCA.csv"),
  list(species = "Maize", file = "data/WorldClim/worldclim_maize_BIO_PCA.csv"),
  list(species = "Sorghum", file = "data/WorldClim/worldclim_sorghum_BIO_PCA.csv")
)

gff_cfg <- list(
  list(species = "Arabidopsis", file = "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/gff3_ref/TAIR10_GFF3_genes.gff"),
  list(species = "Barley", file = "/Users/nirwantandukar/Documents/Research/results/Barley/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.62.gff3"),
  list(species = "Rice", file = "/Users/nirwantandukar/Documents/Research/results/Rice_3001/gff3_ref/osa1_r7.all_models.gff3"),
  list(species = "Maize", file = "/Users/nirwantandukar/Documents/Research/data/GENESPACE/workingDirectory/raw_files/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"),
  list(species = "Sorghum", file = "/Users/nirwantandukar/Documents/Research/data/sorghum_annotation/gene_annotation/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
normalize_chr <- function(x) {
  out <- trimws(as.character(x))
  out <- sub("^9311_", "", out)
  out <- sub("^chr", "", out, ignore.case = TRUE)
  out <- sub("^Chr", "", out)
  out <- sub("^0+", "", out)
  out[out == ""] <- NA_character_
  out
}

chr_rank_df <- function(chr_vec) {
  chr_tbl <- data.frame(chr = unique(chr_vec), stringsAsFactors = FALSE)
  chr_tbl$chr_num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chr_tbl$chr)))
  chr_tbl <- chr_tbl %>% arrange(is.na(chr_num), chr_num, chr)
  chr_tbl$chr_index <- seq_len(nrow(chr_tbl))
  chr_tbl
}

extract_gene_id <- function(attr) {
  id <- sub(".*\\bID=([^;]+).*", "\\1", attr)
  id[id == attr] <- NA_character_
  id <- sub("^gene:", "", id)
  id
}

read_gene_ranges <- function(path) {
  g <- as.data.table(
    read.delim(
      path,
      sep = "\t",
      header = FALSE,
      quote = "",
      comment.char = "#",
      fill = TRUE,
      stringsAsFactors = FALSE
    )
  )
  if (ncol(g) < 9) {
    stop("Invalid GFF/GTF file (fewer than 9 columns): ", path)
  }
  setnames(g, c("chr", "source", "feature", "start", "end", "score", "strand", "phase", "attr"))
  g <- g[feature == "gene"]
  g[, gene := extract_gene_id(attr)]
  g <- g[!is.na(gene) & gene != ""]
  g[, chr := normalize_chr(chr)]
  g <- g[!is.na(chr)]
  g[, `:=`(start = as.integer(start), end = as.integer(end))]
  g <- g[!is.na(start) & !is.na(end)]
  g[, c("source", "feature", "score", "strand", "phase", "attr") := NULL]
  unique(g)
}

annotate_snps_to_genes <- function(snps_df, gene_ranges, species_name, window = 25000L) {
  if (nrow(snps_df) == 0) {
    return(tibble(species = character(), gene = character(), pvalue = numeric()))
  }

  snps <- as.data.table(snps_df)
  snps[, chr := normalize_chr(chr)]
  snps <- snps[!is.na(chr) & !is.na(pos) & !is.na(p) & p > 0]
  if (nrow(snps) == 0) {
    return(tibble(species = character(), gene = character(), pvalue = numeric()))
  }

  snps[, start := pmax(1L, as.integer(pos - window))]
  snps[, end := as.integer(pos + window)]
  snps <- unique(snps, by = c("snp", "chr", "pos", "p", "start", "end"))

  genes <- copy(gene_ranges)
  setkey(genes, chr, start, end)
  setkey(snps, chr, start, end)

  ov <- foverlaps(snps, genes, nomatch = 0L)
  if (nrow(ov) == 0) {
    return(tibble(species = character(), gene = character(), pvalue = numeric()))
  }

  ov[, dist_to_gene := fifelse(
    pos < i.start, i.start - pos,
    fifelse(pos > i.end, pos - i.end, 0L)
  )]

  # Keep closest gene per SNP (ties retained).
  ov <- ov[, .SD[dist_to_gene == min(dist_to_gene)], by = .(snp, chr, pos, p)]

  out <- as_tibble(ov) %>%
    transmute(species = species_name, gene = gene, pvalue = p) %>%
    group_by(species, gene) %>%
    summarise(pvalue = min(pvalue, na.rm = TRUE), .groups = "drop") %>%
    arrange(pvalue, gene)

  out
}

load_gwas <- function(cfg) {
  dat <- read.csv(cfg$file, sep = cfg$sep, check.names = FALSE)

  out <- dat %>%
    transmute(
      species = cfg$species,
      snp = as.character(.data[[cfg$snp_col]]),
      chr = normalize_chr(.data[[cfg$chr_col]]),
      pos = as.numeric(.data[[cfg$pos_col]]),
      p = as.numeric(.data[[cfg$p_col]])
    ) %>%
    filter(!is.na(chr), chr != "", !is.na(pos), !is.na(p), p > 0) %>%
    group_by(species, snp, chr, pos) %>%
    summarise(p = min(p, na.rm = TRUE), .groups = "drop")

  n_total <- nrow(out)
  n_rank <- max(1L, floor(n_total * top_fraction))
  p_cutoff <- out %>% arrange(p) %>% slice(n_rank) %>% pull(p)
  bonf_cutoff <- 0.05 / n_total

  out %>%
    mutate(
      logp = -log10(p),
      p_cutoff = p_cutoff,
      cutoff_logp = -log10(p_cutoff),
      bonf_cutoff = bonf_cutoff,
      bonf_logp = -log10(bonf_cutoff)
    )
}

load_pheno <- function(cfg) {
  dat <- read.csv(cfg$file, check.names = FALSE)
  dat %>%
    transmute(
      species = cfg$species,
      pc1 = as.numeric(PC1)
    ) %>%
    filter(!is.na(pc1))
}

# ------------------------------------------------------------
# Load GWAS + phenotype
# ------------------------------------------------------------
gwas_df <- bind_rows(lapply(gwas_cfg, load_gwas)) %>%
  mutate(species = factor(species, levels = species_order))

pheno_df <- bind_rows(lapply(pheno_cfg, load_pheno)) %>%
  mutate(
    species = factor(species, levels = species_order),
    crop = species_to_crop[as.character(species)]
  )

fill_map <- setNames(crop_colors[unname(species_to_crop[species_order])], species_order)

# ------------------------------------------------------------
# Selection tables
# ------------------------------------------------------------
selection_counts <- gwas_df %>%
  group_by(species) %>%
  summarise(
    n_total = n(),
    n_rank_target = first(max(1L, floor(n() * top_fraction))),
    p_cutoff = first(p_cutoff),
    n_top0p5_snps = sum(p <= first(p_cutoff)),
    bonf_cutoff = first(bonf_cutoff),
    n_bonf_snps = sum(p <= first(bonf_cutoff)),
    .groups = "drop"
  ) %>%
  arrange(species)

top_snps <- gwas_df %>%
  group_by(species) %>%
  filter(p <= first(p_cutoff)) %>%
  ungroup() %>%
  arrange(species, p, chr, pos)

bonf_snps <- gwas_df %>%
  group_by(species) %>%
  filter(p <= first(bonf_cutoff)) %>%
  ungroup() %>%
  arrange(species, p, chr, pos)

write.csv(selection_counts, file.path(res_dir, "pc1_top0p5_and_bonf_selection_counts.csv"), row.names = FALSE)
write.csv(top_snps, file.path(res_dir, "pc1_top0p5_selected_snps.csv"), row.names = FALSE)
write.csv(bonf_snps, file.path(res_dir, "pc1_bonferroni_selected_snps.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Manhattan + phenotype histogram (Fig4)
# ------------------------------------------------------------
prepared <- gwas_df %>%
  group_by(species) %>%
  group_modify(~ {
    x <- .x
    sp <- as.character(.y$species[[1]])
    chr_tbl <- chr_rank_df(x$chr)

    chr_sizes <- x %>%
      group_by(chr) %>%
      summarise(chr_max = max(pos, na.rm = TRUE), .groups = "drop") %>%
      left_join(chr_tbl, by = "chr") %>%
      arrange(chr_index)

    chr_sizes$offset <- c(0, cumsum(head(chr_sizes$chr_max, -1)))

    x2 <- x %>%
      left_join(chr_sizes %>% select(chr, chr_index, offset), by = "chr") %>%
      mutate(
        cum_pos = pos + offset,
        chr_color = ifelse(chr_index %% 2 == 0, "grey65", "black")
      )

    if (identical(sp, "Rice")) {
      x2 <- x2 %>% filter(logp <= 20) %>% mutate(logp_plot = logp)
    } else {
      x2 <- x2 %>% mutate(logp_plot = logp)
    }
    x2
  }) %>%
  ungroup()

cutoffs <- prepared %>%
  group_by(species) %>%
  summarise(p_cutoff = first(p_cutoff), cutoff_logp = first(cutoff_logp), .groups = "drop")

write.csv(cutoffs, file.path(res_dir, "pc1_top0p5_cutoffs.csv"), row.names = FALSE)

p_manhattan <- ggplot(prepared, aes(x = cum_pos, y = logp_plot)) +
  geom_point(aes(color = chr_color), size = 0.12, alpha = 0.65) +
  scale_color_identity() +
  geom_hline(data = cutoffs, aes(yintercept = cutoff_logp), linetype = "dashed", color = "red", linewidth = 0.35) +
  facet_wrap(~ species, ncol = 1, scales = "free") +
  labs(
    title = "Manhattan Plots (WorldClim PC1 GWAS; Top 0.5% SNP Threshold Per Species)",
    x = "Chromosome position",
    y = expression(-log[10](p))
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

p_pheno <- ggplot(pheno_df, aes(x = pc1, fill = species)) +
  geom_histogram(bins = 35, color = "white", linewidth = 0.15) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = fill_map) +
  labs(
    title = "WorldClim PC1 Phenotype Distribution",
    x = "PC1 score",
    y = "Count"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

fig4 <- p_manhattan + p_pheno +
  plot_layout(widths = c(4.5, 1.5)) +
  plot_annotation(tag_levels = "A")

fig4_out <- file.path(fig_dir, "Fig4_manhattan_plus_PC1.png")
ggsave(fig4_out, fig4, width = 18, height = 12, dpi = 300, bg = "white", device = ragg::agg_png)

# ------------------------------------------------------------
# Gene mapping (Bonferroni + top0.5%)
# ------------------------------------------------------------
gene_ranges <- lapply(gff_cfg, function(x) {
  gr <- read_gene_ranges(x$file)
  gr[, species := x$species]
  gr
})
names(gene_ranges) <- vapply(gff_cfg, `[[`, character(1), "species")

top_gene_tbl <- bind_rows(lapply(species_order, function(sp) {
  ss <- top_snps %>%
    filter(species == sp) %>%
    transmute(snp, chr, pos, p)
  annotate_snps_to_genes(ss, gene_ranges[[sp]], sp, window = window_bp)
}))

bonf_gene_tbl <- bind_rows(lapply(species_order, function(sp) {
  ss <- bonf_snps %>%
    filter(species == sp) %>%
    transmute(snp, chr, pos, p)
  annotate_snps_to_genes(ss, gene_ranges[[sp]], sp, window = window_bp)
}))

top_gene_tbl <- top_gene_tbl %>% arrange(species, pvalue, gene)
bonf_gene_tbl <- bonf_gene_tbl %>% arrange(species, pvalue, gene)

write.csv(top_gene_tbl, file.path(res_dir, "pc1_top0p5_gene_table.csv"), row.names = FALSE)
write.csv(bonf_gene_tbl, file.path(res_dir, "pc1_bonferroni_gene_table.csv"), row.names = FALSE)

# Manuscript-ready Bonferroni table requested: species, gene, pvalue
bonf_out <- bonf_gene_tbl %>%
  transmute(species = tolower(species), gene, pvalue)
write.csv(bonf_out, file.path(res_dir, "pc1_bonferroni_top_genes_species_gene_pvalue.csv"), row.names = FALSE)

message("Saved Fig4: ", fig4_out)
message("Saved PC1 GWAS tables in: ", res_dir)
