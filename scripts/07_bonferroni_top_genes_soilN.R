suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# ------------------------------------------------------------
# Bonferroni-significant gene table for soil N GWAS (0-5 cm)
# ------------------------------------------------------------

out_dir <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_cfg <- list(
  list(
    species = "arabidopsis",
    test_file = "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/GWAS_annotate/arabidopsis_TN_0_5_mean.annot_25000bp.tsv",
    test_sep = "tsv",
    test_chr = "CHR",
    test_pos = "BP",
    test_p = "P",
    gene_file = "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/GWAS_annotate/arabidopsis_TN_0_5_mean.annot_25000bp.tsv",
    gene_sep = "tsv",
    gene_col = "closest_gene",
    gene_p = "P"
  ),
  list(
    species = "barley",
    test_file = "/Users/nirwantandukar/Documents/Research/results/Barley/GWAS_annotate/barley_soilN_mod_sub_barley_soilN_gwas_phenotype_gbs_landrace_12129inds_beagle_imputed_isec558kSNP.chrfix.assoc.annot_25000bp.tsv",
    test_sep = "tsv",
    test_chr = "chr",
    test_pos = "ps",
    test_p = "p_wald",
    gene_file = "/Users/nirwantandukar/Documents/Research/results/Barley/GWAS_annotate/barley_soilN_mod_sub_barley_soilN_gwas_phenotype_gbs_landrace_12129inds_beagle_imputed_isec558kSNP.chrfix.assoc.annot_25000bp.tsv",
    gene_sep = "tsv",
    gene_col = "closest_gene",
    gene_p = "p_wald"
  ),
  list(
    species = "maize",
    test_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_maize_LMM.txt",
    test_sep = "tsv",
    test_chr = "chr",
    test_pos = "ps",
    test_p = "p_wald",
    gene_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/annotation_maize_LMM_nitrogen_0_5.csv",
    gene_sep = "csv",
    gene_col = "GeneID",
    gene_p = "PValue"
  ),
  list(
    species = "rice",
    test_file = "/Users/nirwantandukar/Documents/Research/results/Rice_3001/GWAS_annotate/rice_N_mod_sub_rice_gwas_phenotype_TN_rice3000_gwas_qc.snp.assoc.annot_25000bp.tsv",
    test_sep = "tsv",
    test_chr = "chr",
    test_pos = "ps",
    test_p = "p_wald",
    gene_file = "/Users/nirwantandukar/Documents/Research/results/Rice_3001/GWAS_annotate/rice_N_mod_sub_rice_gwas_phenotype_TN_rice3000_gwas_qc.snp.assoc.annot_25000bp.tsv",
    gene_sep = "tsv",
    gene_col = "closest_gene",
    gene_p = "p_wald"
  ),
  list(
    species = "sorghum",
    test_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_sorghum_LMM.txt",
    test_sep = "tsv",
    test_chr = "chr",
    test_pos = "ps",
    test_p = "p_wald",
    gene_file = "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/annotation_sorghum_LMM_nitrogen_0-5_sorghum.csv",
    gene_sep = "csv",
    gene_col = "GeneID",
    gene_p = "PValue"
  )
)

read_any <- function(path, sep_type, cols = NULL) {
  if (!file.exists(path)) stop("Missing file: ", path)
  if (sep_type == "tsv") {
    read_tsv(path, col_types = cols, show_col_types = FALSE, progress = FALSE)
  } else {
    read_csv(path, col_types = cols, show_col_types = FALSE, progress = FALSE)
  }
}

count_tests <- function(cfg) {
  dat <- read_any(cfg$test_file, cfg$test_sep)
  req <- c(cfg$test_chr, cfg$test_pos, cfg$test_p)
  miss <- setdiff(req, names(dat))
  if (length(miss) > 0) stop("Missing required columns in ", cfg$test_file, ": ", paste(miss, collapse = ", "))

  dat %>%
    transmute(
      chr = as.character(.data[[cfg$test_chr]]),
      pos = as.numeric(.data[[cfg$test_pos]]),
      p = as.numeric(.data[[cfg$test_p]])
    ) %>%
    filter(!is.na(chr), str_trim(chr) != "", !is.na(pos), !is.na(p), p > 0) %>%
    distinct(chr, pos) %>%
    nrow()
}

get_sig_genes <- function(cfg, cutoff) {
  dat <- read_any(cfg$gene_file, cfg$gene_sep)
  req <- c(cfg$gene_col, cfg$gene_p)
  miss <- setdiff(req, names(dat))
  if (length(miss) > 0) stop("Missing required columns in ", cfg$gene_file, ": ", paste(miss, collapse = ", "))

  out <- dat %>%
    transmute(
      species = cfg$species,
      gene = as.character(.data[[cfg$gene_col]]),
      pvalue = as.numeric(.data[[cfg$gene_p]])
    ) %>%
    filter(
      !is.na(gene),
      str_trim(gene) != "",
      toupper(str_trim(gene)) != "NA",
      is.finite(pvalue),
      pvalue > 0,
      pvalue <= cutoff
    )

  if (nrow(out) == 0) {
    return(tibble(species = character(), gene = character(), pvalue = numeric()))
  }

  out %>%
    group_by(species, gene) %>%
    summarise(pvalue = min(pvalue), .groups = "drop")
}

summary_rows <- vector("list", length(species_cfg))
sig_rows <- vector("list", length(species_cfg))

for (i in seq_along(species_cfg)) {
  cfg <- species_cfg[[i]]
  n_tests <- count_tests(cfg)
  cutoff <- 0.05 / n_tests
  sig <- get_sig_genes(cfg, cutoff)

  summary_rows[[i]] <- tibble(
    species = cfg$species,
    n_tests = n_tests,
    bonferroni_cutoff = cutoff,
    n_significant_genes = nrow(sig)
  )
  sig_rows[[i]] <- sig
}

summary_tbl <- bind_rows(summary_rows) %>% arrange(species)
sig_tbl <- bind_rows(sig_rows) %>% arrange(species, pvalue, gene)

summary_out <- file.path(out_dir, "soilN_bonferroni_gene_summary.csv")
sig_out <- file.path(out_dir, "soilN_bonferroni_significant_genes_table.csv")

write_csv(summary_tbl, summary_out)
write_csv(sig_tbl, sig_out)

message("Saved: ", summary_out)
message("Saved: ", sig_out)
message("Total significant genes: ", nrow(sig_tbl))
