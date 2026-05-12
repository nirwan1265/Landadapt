suppressPackageStartupMessages({
  library(SNPRelate)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[length(idx)] == length(args)) stop("Missing value for ", flag)
  args[idx[length(idx)] + 1]
}

vcf_file <- get_arg("--vcf", "/Users/nirwantandukar/Documents/Research/data/Rice/rice3000_main12_chrnum.vcf.gz")
chr_keep <- get_arg("--chr", "1")
out_dir <- get_arg("--outdir", file.path("results", "rice_snprelate_chr1"))
max_random_snps <- as.integer(get_arg("--max-random-snps", "100000"))
ld_threshold <- as.numeric(get_arg("--ld-threshold", "0.2"))
maf_cutoff <- as.numeric(get_arg("--maf", "0.05"))
missing_rate <- as.numeric(get_arg("--missing-rate", "0.1"))
threads <- as.integer(get_arg("--threads", "4"))
seed <- as.integer(get_arg("--seed", "2026"))

if (!file.exists(vcf_file)) stop("VCF not found: ", vcf_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)

chr_vcf <- file.path(out_dir, paste0("chr", chr_keep, ".vcf.gz"))
gds_file <- file.path(out_dir, paste0("chr", chr_keep, ".gds"))
pcs_csv <- file.path(out_dir, paste0("rice_chr", chr_keep, "_genotype_PCs.csv"))
var_csv <- file.path(out_dir, paste0("rice_chr", chr_keep, "_genotype_PCs_variance.csv"))
plot_png <- file.path(out_dir, paste0("rice_chr", chr_keep, "_PCA_PC1_PC2.png"))

message("Step 1/5: Subset VCF to chromosome ", chr_keep)
if (!file.exists(chr_vcf)) {
  cmd_subset <- sprintf(
    "bcftools view -r %s -Oz -o %s %s",
    shQuote(chr_keep), shQuote(chr_vcf), shQuote(vcf_file)
  )
  status <- system(cmd_subset, ignore.stdout = TRUE, ignore.stderr = FALSE)
  if (!identical(status, 0L)) stop("bcftools subset failed")
}
cmd_index <- sprintf("bcftools index -f -c %s", shQuote(chr_vcf))
status <- system(cmd_index, ignore.stdout = TRUE, ignore.stderr = FALSE)
if (!identical(status, 0L)) stop("bcftools index failed")

message("Step 2/5: Convert chr VCF -> GDS")
snpgdsVCF2GDS(
  vcf.fn = chr_vcf,
  out.fn = gds_file,
  method = "biallelic.only",
  verbose = TRUE
)

message("Step 3/5: Random SNP subset")
genofile <- snpgdsOpen(gds_file, readonly = FALSE)
on.exit(try(snpgdsClose(genofile), silent = TRUE), add = TRUE)

snp_ids <- read.gdsn(index.gdsn(genofile, "snp.id"))
n_total <- length(snp_ids)
if (n_total == 0) stop("No SNPs found on chromosome ", chr_keep)

if (n_total > max_random_snps) {
  snp_subset <- sample(snp_ids, max_random_snps)
} else {
  snp_subset <- snp_ids
}

message("Total SNPs on chr", chr_keep, ": ", format(n_total, big.mark = ","))
message("Random SNPs used before LD pruning: ", format(length(snp_subset), big.mark = ","))

message("Step 4/5: LD pruning")
set.seed(seed)
snpset <- snpgdsLDpruning(
  genofile,
  snp.id = snp_subset,
  autosome.only = FALSE,
  method = "corr",
  slide.max.n = 500L,
  ld.threshold = ld_threshold,
  maf = maf_cutoff,
  missing.rate = missing_rate,
  num.thread = threads,
  verbose = TRUE
)
pruned_snps <- unlist(snpset, use.names = FALSE)
if (length(pruned_snps) == 0) stop("No SNPs left after LD pruning; try lowering maf or increasing missing-rate")
message("SNPs after LD pruning: ", format(length(pruned_snps), big.mark = ","))

message("Step 5/5: PCA")
pca <- snpgdsPCA(
  genofile,
  snp.id = pruned_snps,
  num.thread = threads,
  autosome.only = FALSE,
  verbose = TRUE
)

pc_df <- tibble(
  sample_id = pca$sample.id,
  PC1 = pca$eigenvect[, 1],
  PC2 = pca$eigenvect[, 2],
  PC3 = pca$eigenvect[, 3]
)

var_df <- tibble(
  PC = paste0("PC", seq_along(pca$varprop)),
  variance_explained = pca$varprop
)

write_csv(pc_df, pcs_csv)
write_csv(var_df, var_csv)

p <- ggplot(pc_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 1.2, alpha = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Rice Genotype PCA (Chr ", chr_keep, ")"),
    subtitle = paste0(
      "Random SNPs before LD pruning: ", format(length(snp_subset), big.mark = ","), "; after pruning: ",
      format(length(pruned_snps), big.mark = ",")
    ),
    x = paste0("PC1 (", round(100 * pca$varprop[1], 2), "%)"),
    y = paste0("PC2 (", round(100 * pca$varprop[2], 2), "%)")
  )
ggsave(plot_png, p, width = 7, height = 5, dpi = 300, bg = "white", device = ragg::agg_png)

message("Done.")
message("PC file: ", pcs_csv)
message("Variance file: ", var_csv)
message("Plot: ", plot_png)
