suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------
# Config
# -----------------------------
fig2_dir <- "Figs/Fig2"
dir.create(fig2_dir, recursive = TRUE, showWarnings = FALSE)

top_fraction <- 0.005

crop_colors <- c(
  maize = "#E69F00",
  sorghum = "#D55E00",
  rice = "#CC79A7",
  arabidopsis = "#009E73",
  barley = "#0072B2",
  pearl_millet = "#F0E442"
)

species_order <- c("Arabidopsis", "Barley", "Rice", "Maize", "Sorghum")

gwas_cfg <- list(
  list(species = "Arabidopsis", file = "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/GWAS_annotate/arabidopsis_TN_0_5_mean.annot_25000bp.tsv", sep = "\t", chr_col = "CHR", pos_col = "BP", p_col = "P"),
  list(species = "Barley", file = "/Users/nirwantandukar/Documents/Research/data/Barley/GWAS_results/barley_soilN_mod_sub_barley_soilN_gwas_phenotype_gbs_landrace_12129inds_beagle_imputed_isec558kSNP.assoc.txt", sep = "\t", chr_col = "chr", pos_col = "ps", p_col = "p_wald"),
  list(species = "Rice", file = "/Users/nirwantandukar/Documents/Research/results/Rice_3001/GWAS_annotate/rice_N_mod_sub_rice_gwas_phenotype_TN_rice3000_gwas_qc.snp.assoc.annot_25000bp.tsv", sep = "\t", chr_col = "chr", pos_col = "ps", p_col = "p_wald"),
  list(species = "Maize", file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_maize_LMM.txt", sep = "\t", chr_col = "chr", pos_col = "ps", p_col = "p_wald"),
  list(species = "Sorghum", file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_sorghum_LMM.txt", sep = "\t", chr_col = "chr", pos_col = "ps", p_col = "p_wald")
)

pheno_cfg <- list(
  list(species = "Arabidopsis", file = "data/TN/arabidopsis_N_values.csv"),
  list(species = "Barley", file = "data/TN/barley_N_values.csv"),
  list(species = "Rice", file = "data/TN/rice_N_values.csv"),
  list(species = "Maize", file = "data/TN/maize_N_values.csv"),
  list(species = "Sorghum", file = "data/TN/sorghum_N_values.csv")
)

species_to_crop <- c(
  Arabidopsis = "arabidopsis",
  Barley = "barley",
  Rice = "rice",
  Maize = "maize",
  Sorghum = "sorghum"
)

# -----------------------------
# Helpers
# -----------------------------
normalize_chr <- function(x) {
  out <- trimws(as.character(x))
  out <- sub("^9311_", "", out)
  out <- sub("^chr", "", out, ignore.case = TRUE)
  out <- sub("^Chr", "", out)
  out
}

chr_rank_df <- function(chr_vec) {
  chr_tbl <- data.frame(chr = unique(chr_vec), stringsAsFactors = FALSE)
  chr_tbl$chr_num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chr_tbl$chr)))
  chr_tbl <- chr_tbl %>% arrange(is.na(chr_num), chr_num, chr)
  chr_tbl$chr_index <- seq_len(nrow(chr_tbl))
  chr_tbl
}

load_gwas <- function(cfg) {
  dat <- read.csv(cfg$file, sep = cfg$sep, check.names = FALSE)
  out <- dat %>%
    transmute(
      species = cfg$species,
      chr = normalize_chr(.data[[cfg$chr_col]]),
      pos = as.numeric(.data[[cfg$pos_col]]),
      p = as.numeric(.data[[cfg$p_col]])
    ) %>%
    filter(!is.na(chr), chr != "", !is.na(pos), !is.na(p), p > 0) %>%
    group_by(species, chr, pos) %>%
    summarise(p = min(p, na.rm = TRUE), .groups = "drop")

  n_total <- nrow(out)
  n_rank <- max(1L, floor(n_total * top_fraction))
  p_cutoff <- out %>% arrange(p) %>% slice(n_rank) %>% pull(p)

  out <- out %>%
    mutate(
      logp = -log10(p),
      p_cutoff = p_cutoff,
      cutoff_logp = -log10(p_cutoff)
    )

  out
}

load_pheno <- function(cfg) {
  dat <- read.csv(cfg$file, check.names = FALSE)
  out <- dat %>%
    mutate(
      species = cfg$species,
      log_n = ifelse(!is.na(log_n), log_n, ifelse(n_value > 0, log10(n_value), NA_real_))
    ) %>%
    filter(!is.na(log_n)) %>%
    select(species, log_n)
  out
}

# -----------------------------
# Prepare GWAS data
# -----------------------------
gwas_df <- bind_rows(lapply(gwas_cfg, load_gwas)) %>%
  mutate(species = factor(species, levels = species_order))

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
      # Drop extreme rice outliers entirely for cleaner visualization.
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

write.csv(cutoffs, file.path(fig2_dir, "fig2_top0p5pct_cutoffs.csv"), row.names = FALSE)

# -----------------------------
# Prepare phenotype data
# -----------------------------
pheno_df <- bind_rows(lapply(pheno_cfg, load_pheno)) %>%
  mutate(
    species = factor(species, levels = species_order),
    crop = species_to_crop[as.character(species)]
  )

fill_map <- setNames(crop_colors[unname(species_to_crop[species_order])], species_order)

# -----------------------------
# Plot panels
# -----------------------------
p_manhattan <- ggplot(prepared, aes(x = cum_pos, y = logp_plot)) +
  geom_point(aes(color = chr_color), size = 0.12, alpha = 0.65) +
  scale_color_identity() +
  geom_hline(data = cutoffs, aes(yintercept = cutoff_logp), linetype = "dashed", color = "red", linewidth = 0.35) +
  facet_wrap(~ species, ncol = 1, scales = "free") +
  labs(
    title = "Manhattan Plots (Top 0.5% SNP Threshold Per Species)",
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

p_pheno <- ggplot(pheno_df, aes(x = log_n, fill = species)) +
  geom_histogram(bins = 35, color = "white", linewidth = 0.15) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = fill_map) +
  labs(
    title = "Nitrogen Phenotype Distribution",
    x = expression(log[10](N)),
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

fig2 <- p_manhattan + p_pheno +
  plot_layout(widths = c(4.5, 1.5)) +
  plot_annotation(tag_levels = "A")

fig2_out <- file.path(fig2_dir, "Fig2_manhattan_plus_nitrogen.png")
ggsave(fig2_out, fig2, width = 18, height = 12, dpi = 300, bg = "white", device = ragg::agg_png)

message("Saved Fig2: ", fig2_out)
