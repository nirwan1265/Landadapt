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
args <- commandArgs(trailingOnly = TRUE)
trait_input <- if (length(args) >= 1) trimws(args[[1]]) else "PC1"

normalize_trait <- function(x) {
  key <- toupper(gsub("[^A-Za-z0-9]+", "", x))
  if (key %in% c("PC1", "WORLDCLIMPC1")) return("PC1")
  if (key %in% c("PC2", "WORLDCLIMPC2")) return("PC2")
  if (key %in% c("PC3", "WORLDCLIMPC3")) return("PC3")
  if (key %in% c("CEC")) return("CEC")
  if (key %in% c("PH", "SOILPH")) return("PH")
  if (key %in% c("SOILN", "TN", "N")) return("SOILN")
  if (key %in% c("AMROOTSCOLONIZED", "AMROOTS", "AMROOT")) return("AM_ROOTS_COLONIZED")
  if (key %in% c("AMRELABUNDANCECOLONIZATION", "AMRELABUNDANCECOLONIZATION", "AMREL", "AMABUNDANCE")) return("AM_REL_ABUNDANCE_COLONIZATION")
  if (key %in% c("ARIDITY", "ARIDITYINDEX", "AI")) return("ARIDITY_INDEX")
  NA_character_
}

trait <- normalize_trait(trait_input)
if (is.na(trait)) {
  stop(
    "Unsupported trait: ", trait_input,
    "\nSupported: PC1, PC2, PC3, CEC, PH, SOILN, AM_ROOTS_COLONIZED, AM_REL_ABUNDANCE_COLONIZATION, ARIDITY_INDEX"
  )
}

trait_cfg <- list(
  PC1 = list(
    gwas_token = "PC1",
    result_dir = "results/pc1_worldclim",
    out_prefix = "pc1_worldclim",
    fig_suffix = "PC1",
    pheno_file = function(sp) file.path("data/WorldClim", paste0("worldclim_", sp, "_BIO_PCA.csv")),
    pheno_col = "PC1",
    manhattan_title = "Manhattan Plots (WorldClim PC1 GWAS; Top 0.5% SNP Threshold Per Species)",
    pheno_title = "WorldClim PC1 Phenotype Distribution",
    pheno_x = "PC1 score"
  ),
  PC2 = list(
    gwas_token = "PC2",
    result_dir = "results/pc2_worldclim",
    out_prefix = "pc2_worldclim",
    fig_suffix = "PC2",
    pheno_file = function(sp) file.path("data/WorldClim", paste0("worldclim_", sp, "_BIO_PCA.csv")),
    pheno_col = "PC2",
    manhattan_title = "Manhattan Plots (WorldClim PC2 GWAS; Top 0.5% SNP Threshold Per Species)",
    pheno_title = "WorldClim PC2 Phenotype Distribution",
    pheno_x = "PC2 score"
  ),
  PC3 = list(
    gwas_token = "PC3",
    result_dir = "results/pc3_worldclim",
    out_prefix = "pc3_worldclim",
    fig_suffix = "PC3",
    pheno_file = function(sp) file.path("data/WorldClim", paste0("worldclim_", sp, "_BIO_PCA.csv")),
    pheno_col = "PC3",
    manhattan_title = "Manhattan Plots (WorldClim PC3 GWAS; Top 0.5% SNP Threshold Per Species)",
    pheno_title = "WorldClim PC3 Phenotype Distribution",
    pheno_x = "PC3 score"
  ),
  CEC = list(
    gwas_token = "cec",
    result_dir = "results/cec",
    out_prefix = "cec",
    fig_suffix = "CEC",
    pheno_file = function(sp) file.path("data/Cec", paste0(sp, "_Cec_values.csv")),
    pheno_col = "cec_value",
    manhattan_title = "Manhattan Plots (CEC GWAS; Top 0.5% SNP Threshold Per Species)",
    pheno_title = "Soil CEC Phenotype Distribution",
    pheno_x = "CEC value"
  ),
  PH = list(
    gwas_token = "ph",
    result_dir = "results/ph",
    out_prefix = "ph",
    fig_suffix = "pH",
    pheno_file = function(sp) file.path("data/pH", paste0(sp, "_pH_values.csv")),
    pheno_col = "ph_value",
    manhattan_title = "Manhattan Plots (pH GWAS; Top 0.5% SNP Threshold Per Species)",
    pheno_title = "Soil pH Phenotype Distribution",
    pheno_x = "pH value"
  ),
  SOILN = list(
    gwas_token = "soilN",
    result_dir = "results/soilN",
    out_prefix = "soilN",
    fig_suffix = "soilN",
    pheno_file = function(sp) file.path("data/TN", paste0(sp, "_N_values.csv")),
    pheno_col = "n_value",
    manhattan_title = "Manhattan Plots (Soil N GWAS; Top 0.5% SNP Threshold Per Species)",
    pheno_title = "Soil Nitrogen (0-5 cm) Phenotype Distribution",
    pheno_x = "Soil N value"
  ),
  AM_ROOTS_COLONIZED = list(
    gwas_token = "am_roots_colonized",
    result_dir = "results/am_roots_colonized",
    out_prefix = "am_roots_colonized",
    fig_suffix = "AM_roots_colonized",
    pheno_file = function(sp) file.path("data/AM_roots_colonized", paste0(sp, "_AM_roots_colonized_values.csv")),
    pheno_col = "am_roots_colonized_value",
    manhattan_title = "Manhattan Plots (AM Roots Colonized GWAS; Top 0.5% SNP Threshold Per Species)",
    pheno_title = "AM Roots Colonized Phenotype Distribution",
    pheno_x = "AM roots colonized"
  ),
  AM_REL_ABUNDANCE_COLONIZATION = list(
    gwas_token = "am_rel_abundance_colonization",
    result_dir = "results/am_rel_abundance_colonization",
    out_prefix = "am_rel_abundance_colonization",
    fig_suffix = "AM_rel_abundance_colonization",
    pheno_file = function(sp) file.path("data/AM_rel_abundance_colonization", paste0(sp, "_AM_rel_abundance_colonization_values.csv")),
    pheno_col = "am_rel_abundance_colonization_value",
    manhattan_title = "Manhattan Plots (AM Relative Abundance GWAS; Top 0.5% SNP Threshold Per Species)",
    pheno_title = "AM Relative Abundance Phenotype Distribution",
    pheno_x = "AM relative abundance"
  ),
  ARIDITY_INDEX = list(
    gwas_token = "aridity_index",
    result_dir = "results/aridity_index",
    out_prefix = "aridity_index",
    fig_suffix = "aridity_index",
    pheno_file = function(sp) file.path("data/Aridity", paste0(sp, "_aridity_values.csv")),
    pheno_col = "aridity_index",
    manhattan_title = "Manhattan Plots (Aridity Index GWAS; Top 0.5% SNP Threshold Per Species)",
    pheno_title = "Aridity Index Phenotype Distribution",
    pheno_x = "Aridity index"
  )
)

cfg <- trait_cfg[[trait]]

fig_dir <- "Figs/Fig4"
res_dir <- cfg$result_dir
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

gwas_base <- "/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt"
make_gwas_cfg <- function(sp) {
  sp_lc <- species_to_crop[[sp]]
  prefix <- if (sp_lc == "arabidopsis") "AT" else sp_lc
  list(
    species = sp,
    file = file.path(gwas_base, sp_lc, paste0(prefix, "_", cfg$gwas_token, ".txt")),
    sep = "\t",
    chr_col = "chr",
    pos_col = "ps",
    p_col = "p_wald",
    snp_col = "rs"
  )
}
gwas_cfg <- lapply(species_order, make_gwas_cfg)

make_pheno_cfg <- function(sp) {
  sp_lc <- species_to_crop[[sp]]
  list(species = sp, file = cfg$pheno_file(sp_lc), col = cfg$pheno_col)
}
pheno_cfg <- lapply(species_order, make_pheno_cfg)

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

empty_gwas <- function(sp) {
  tibble(
    species = character(),
    snp = character(),
    chr = character(),
    pos = numeric(),
    p = numeric(),
    logp = numeric(),
    p_cutoff = numeric(),
    cutoff_logp = numeric(),
    bonf_cutoff = numeric(),
    bonf_logp = numeric()
  )
}

load_gwas <- function(cfg) {
  if (!file.exists(cfg$file)) {
    message("[WARN] Missing GWAS file for ", cfg$species, ": ", cfg$file)
    return(empty_gwas(cfg$species))
  }

  dat <- read.csv(cfg$file, sep = cfg$sep, check.names = FALSE)
  resolve_col <- function(nms, candidates) {
    idx <- match(tolower(candidates), tolower(nms))
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) return(NA_character_)
    nms[idx[[1]]]
  }

  snp_col <- resolve_col(names(dat), c(cfg$snp_col, "SNP", "snp", "Marker", "marker", "rsid", "RS"))
  chr_col <- resolve_col(names(dat), c(cfg$chr_col, "Chr", "CHR", "CHROM", "chrom", "chromosome"))
  pos_col <- resolve_col(names(dat), c(cfg$pos_col, "Pos", "POS", "position", "bp", "BP"))
  p_col <- resolve_col(names(dat), c(cfg$p_col, "P.value", "p.value", "P", "pvalue", "p_value", "p_lrt", "p_score"))

  if (any(is.na(c(snp_col, chr_col, pos_col, p_col)))) {
    message(
      "[WARN] Could not resolve GWAS columns in ", cfg$file,
      " (snp=", snp_col, ", chr=", chr_col, ", pos=", pos_col, ", p=", p_col, ")"
    )
    return(empty_gwas(cfg$species))
  }

  out <- dat %>%
    transmute(
      species = cfg$species,
      snp = as.character(.data[[snp_col]]),
      chr = normalize_chr(.data[[chr_col]]),
      pos = as.numeric(.data[[pos_col]]),
      p = as.numeric(.data[[p_col]])
    ) %>%
    filter(!is.na(chr), chr != "", !is.na(pos), !is.na(p), p > 0) %>%
    group_by(species, snp, chr, pos) %>%
    summarise(p = min(p, na.rm = TRUE), .groups = "drop")

  n_total <- nrow(out)
  if (n_total == 0) {
    message("[WARN] No usable SNP rows after filtering in ", cfg$file)
    return(empty_gwas(cfg$species))
  }

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
  if (!file.exists(cfg$file)) {
    message("[WARN] Missing phenotype file for ", cfg$species, ": ", cfg$file)
    return(tibble(species = character(), trait_value = numeric()))
  }

  dat <- read.csv(cfg$file, check.names = FALSE)
  if (!cfg$col %in% names(dat)) {
    message("[WARN] Missing phenotype column '", cfg$col, "' in ", cfg$file)
    return(tibble(species = character(), trait_value = numeric()))
  }

  dat %>%
    transmute(
      species = cfg$species,
      trait_value = as.numeric(.data[[cfg$col]])
    ) %>%
    filter(!is.na(trait_value))
}

# ------------------------------------------------------------
# Load GWAS + phenotype
# ------------------------------------------------------------
gwas_df <- bind_rows(lapply(gwas_cfg, load_gwas))
if (nrow(gwas_df) == 0) {
  stop("No GWAS rows available for trait ", trait, ". Check input files.")
}

gwas_df <- gwas_df %>% mutate(species = factor(species, levels = species_order))

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
    n_rank_target = max(1L, floor(n_total * top_fraction)),
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

write.csv(selection_counts, file.path(res_dir, paste0(cfg$out_prefix, "_top0p5_and_bonf_selection_counts.csv")), row.names = FALSE)
write.csv(top_snps, file.path(res_dir, paste0(cfg$out_prefix, "_top0p5_selected_snps.csv")), row.names = FALSE)
write.csv(bonf_snps, file.path(res_dir, paste0(cfg$out_prefix, "_bonferroni_selected_snps.csv")), row.names = FALSE)

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

write.csv(cutoffs, file.path(res_dir, paste0(cfg$out_prefix, "_top0p5_cutoffs.csv")), row.names = FALSE)

p_manhattan <- ggplot(prepared, aes(x = cum_pos, y = logp_plot)) +
  geom_point(aes(color = chr_color), size = 0.12, alpha = 0.65) +
  scale_color_identity() +
  geom_hline(data = cutoffs, aes(yintercept = cutoff_logp), linetype = "dashed", color = "red", linewidth = 0.35) +
  facet_wrap(~ species, ncol = 1, scales = "free") +
  labs(
    title = cfg$manhattan_title,
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

p_pheno <- ggplot(pheno_df, aes(x = trait_value, fill = species)) +
  geom_histogram(bins = 35, color = "white", linewidth = 0.15) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = fill_map) +
  labs(
    title = cfg$pheno_title,
    x = cfg$pheno_x,
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

fig4_out <- file.path(fig_dir, paste0("Fig4_manhattan_plus_", cfg$fig_suffix, ".png"))
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

write.csv(top_gene_tbl, file.path(res_dir, paste0(cfg$out_prefix, "_top0p5_gene_table.csv")), row.names = FALSE)
write.csv(bonf_gene_tbl, file.path(res_dir, paste0(cfg$out_prefix, "_bonferroni_gene_table.csv")), row.names = FALSE)

# Manuscript-ready Bonferroni table requested: species, gene, pvalue
bonf_out <- bonf_gene_tbl %>%
  transmute(species = tolower(species), gene, pvalue)
write.csv(bonf_out, file.path(res_dir, paste0(cfg$out_prefix, "_bonferroni_top_genes_species_gene_pvalue.csv")), row.names = FALSE)

# Also save top0.5% in species,gene,pvalue format for orthogroup mapping.
top_out <- top_gene_tbl %>%
  transmute(species = tolower(species), gene, pvalue)
write.csv(top_out, file.path(res_dir, paste0(cfg$out_prefix, "_top0p5_genes_species_gene_pvalue.csv")), row.names = FALSE)

message("Saved Fig4: ", fig4_out)
message("Saved trait GWAS tables in: ", res_dir)
