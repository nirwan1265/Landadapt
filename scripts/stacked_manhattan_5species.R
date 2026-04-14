suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------
out_dir <- "data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

top_fraction <- 0.005

species_cfg <- list(
  list(
    species = "Arabidopsis",
    file = "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/GWAS_annotate/arabidopsis_TN_0_5_mean.annot_25000bp.tsv",
    sep = "\t",
    chr_col = "CHR",
    pos_col = "BP",
    p_col = "P"
  ),
  list(
    species = "Rice",
    file = "/Users/nirwantandukar/Documents/Research/results/Rice_3001/GWAS_annotate/rice_N_mod_sub_rice_gwas_phenotype_TN_rice3000_gwas_qc.snp.assoc.annot_25000bp.tsv",
    sep = "\t",
    chr_col = "chr",
    pos_col = "ps",
    p_col = "p_wald"
  ),
  list(
    species = "Maize",
    file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_maize_LMM.txt",
    sep = "\t",
    chr_col = "chr",
    pos_col = "ps",
    p_col = "p_wald"
  ),
  list(
    species = "Sorghum",
    file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_sorghum_LMM.txt",
    sep = "\t",
    chr_col = "chr",
    pos_col = "ps",
    p_col = "p_wald"
  ),
  list(
    species = "Barley",
    file = "/Users/nirwantandukar/Documents/Research/data/Barley/GWAS_results/barley_soilN_mod_sub_barley_soilN_gwas_phenotype_gbs_landrace_12129inds_beagle_imputed_isec558kSNP.assoc.txt",
    sep = "\t",
    chr_col = "chr",
    pos_col = "ps",
    p_col = "p_wald"
  )
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
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

load_species <- function(cfg) {
  if (!file.exists(cfg$file)) {
    warning("Missing file for ", cfg$species, ": ", cfg$file)
    return(NULL)
  }

  dat <- read.csv(cfg$file, sep = cfg$sep, check.names = FALSE)

  req <- c(cfg$chr_col, cfg$pos_col, cfg$p_col)
  miss <- setdiff(req, names(dat))
  if (length(miss) > 0) {
    stop("Missing required columns in ", cfg$file, ": ", paste(miss, collapse = ", "))
  }

  out <- dat %>%
    transmute(
      species = cfg$species,
      chr = normalize_chr(.data[[cfg$chr_col]]),
      pos = as.numeric(.data[[cfg$pos_col]]),
      p = as.numeric(.data[[cfg$p_col]])
    ) %>%
    filter(!is.na(chr), chr != "", !is.na(pos), !is.na(p), p > 0)

  # Some annotation tables repeat SNP rows per nearby gene.
  # Keep one point per chr-position with minimum p-value.
  out <- out %>%
    group_by(species, chr, pos) %>%
    summarise(p = min(p, na.rm = TRUE), .groups = "drop")

  # Rank cutoff from total SNPs: p at rank floor(top_fraction * n_total)
  n_total <- nrow(out)
  n_rank <- if (n_total > 0) max(1L, floor(n_total * top_fraction)) else 0L
  p_cutoff <- if (n_rank > 0) out %>% arrange(p) %>% slice(n_rank) %>% pull(p) else NA_real_
  selected <- if (is.finite(p_cutoff)) out %>% filter(p <= p_cutoff) else out[0, ]

  attr(out, "selection_stats") <- data.frame(
    species = cfg$species,
    n_total = n_total,
    n_rank_target = n_rank,
    n_selected = nrow(selected),
    p_cutoff = p_cutoff,
    p_min = if (nrow(selected) > 0) min(selected$p, na.rm = TRUE) else NA_real_,
    p_max = if (nrow(selected) > 0) max(selected$p, na.rm = TRUE) else NA_real_
  )
  attr(out, "selected_snps") <- selected

  out
}

# ------------------------------------------------------------
# Load and prepare all species
# ------------------------------------------------------------
all_list <- lapply(species_cfg, load_species)
all_df <- bind_rows(all_list)
stats_df <- bind_rows(lapply(all_list, function(x) attr(x, "selection_stats")))

if (nrow(all_df) == 0) stop("No GWAS points loaded.")
message("Total points loaded (full Manhattan): ", nrow(all_df))
if (nrow(stats_df) > 0) {
  message("Selection summary (top 0.5% rank cutoff from total SNPs):")
  print(stats_df)
}

stats_csv <- file.path(out_dir, "manhattan_5species_selection_counts_top0p5pct_rankcutoff.csv")
write.csv(stats_df, stats_csv, row.names = FALSE)

selected_df <- bind_rows(lapply(all_list, function(x) attr(x, "selected_snps")))
selected_csv <- file.path(out_dir, "manhattan_5species_selected_snps_top0p5pct_rankcutoff.csv")
write.csv(selected_df %>% arrange(species, p), selected_csv, row.names = FALSE)

prepared <- all_df %>%
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
        logp = -log10(p),
        chr_color = ifelse(chr_index %% 2 == 0, "grey65", "black")
      )
    if (identical(sp, "Rice")) {
      x2 <- x2 %>% mutate(logp_plot = pmin(logp, 20))
    } else {
      x2 <- x2 %>% mutate(logp_plot = logp)
    }

    x2
  }) %>%
  ungroup()

hline_df <- stats_df %>%
  mutate(cutoff_logp = -log10(p_cutoff))

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(prepared, aes(x = cum_pos, y = logp_plot)) +
  geom_point(aes(color = chr_color), size = 0.15, alpha = 0.65) +
  scale_color_identity() +
  geom_hline(
    data = hline_df,
    aes(yintercept = cutoff_logp),
    linetype = "dashed",
    color = "red",
    linewidth = 0.35
  ) +
  facet_wrap(~ species, ncol = 1, scales = "free") +
  labs(
    title = "Stacked Manhattan Plots (5 Species)",
    subtitle = "Full Manhattan (all SNPs), red dashed line = top 0.5% SNP rank p-value cutoff",
    x = "Chromosome position",
    y = expression(-log[10](p))
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

out_png <- file.path(out_dir, "manhattan_5species_stacked_full_top0p5pct_threshold.png")

ggsave(out_png, p, width = 14, height = 12, dpi = 300, bg = "white")

message("Saved: ", out_png)
message("Saved: ", stats_csv)
message("Saved: ", selected_csv)
