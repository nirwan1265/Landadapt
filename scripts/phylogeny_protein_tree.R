#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(msa)
  library(Biostrings)
  library(phangorn)
  library(ape)
})

parse_args <- function(argv) {
  defaults <- list(
    input_dir = "data/fasta_protein",
    out_dir = "results/phylogeny",
    bootstrap = 100L,
    seed = 2026L
  )

  if (length(argv) == 0) {
    return(defaults)
  }

  for (arg in argv) {
    if (!grepl("^--[^=]+=.*$", arg)) {
      stop("Arguments must be in --key=value format. Got: ", arg, call. = FALSE)
    }
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- kv[2]
    if (key == "input_dir") defaults$input_dir <- val
    if (key == "out_dir") defaults$out_dir <- val
    if (key == "bootstrap") defaults$bootstrap <- as.integer(val)
    if (key == "seed") defaults$seed <- as.integer(val)
  }

  defaults
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

fasta_files <- sort(list.files(
  opts$input_dir,
  full.names = TRUE,
  pattern = "\\.(txt|fa|fasta|faa)$"
))
if (length(fasta_files) < 3) {
  stop("Need at least 3 FASTA files to infer a phylogeny.", call. = FALSE)
}

seqs <- do.call(c, lapply(fasta_files, readAAStringSet))
if (length(seqs) < 3) {
  stop("Need at least 3 sequences to infer a phylogeny.", call. = FALSE)
}

cat("Loaded", length(seqs), "protein sequences\n")

alignment <- msa(seqs, method = "Muscle", type = "protein")
aln_set <- as(alignment, "AAStringSet")
alignment_path <- file.path(opts$out_dir, "orthogroup_protein_alignment.fasta")
writeXStringSet(aln_set, filepath = alignment_path, format = "fasta")
cat("Alignment length:", unique(width(aln_set)), "aa\n")

phyd <- phyDat(as.matrix(aln_set), type = "AA")
dist_matrix <- dist.ml(phyd, model = "JTT")
tree_nj <- NJ(dist_matrix)

fit <- pml(tree_nj, data = phyd)
fit_opt <- optim.pml(
  fit,
  model = "JTT",
  optInv = TRUE,
  optGamma = TRUE,
  rearrangement = "stochastic",
  control = pml.control(trace = 0)
)

set.seed(opts$seed)
bs <- bootstrap.pml(
  fit_opt,
  bs = opts$bootstrap,
  optNni = TRUE,
  multicore = FALSE,
  control = pml.control(trace = 0)
)
tree_bs <- plotBS(fit_opt$tree, bs, p = 50, type = "none")

ml_tree_path <- file.path(opts$out_dir, "orthogroup_protein_tree_ml.newick")
ml_bs_tree_path <- file.path(opts$out_dir, "orthogroup_protein_tree_ml_bootstrap.newick")
pdf_path <- file.path(opts$out_dir, "orthogroup_protein_tree_ml_bootstrap.pdf")
summary_path <- file.path(opts$out_dir, "run_summary.txt")

write.tree(fit_opt$tree, file = ml_tree_path)
write.tree(tree_bs, file = ml_bs_tree_path)

pdf(pdf_path, width = 9, height = 6)
plotBS(
  fit_opt$tree,
  bs,
  p = 50,
  cex = 0.9,
  main = sprintf("Protein phylogeny (JTT ML, %d bootstrap)", opts$bootstrap)
)
add.scale.bar()
dev.off()

writeLines(
  c(
    paste0("n_sequences: ", length(seqs)),
    paste0("sequence_names: ", paste(names(seqs), collapse = ", ")),
    paste0("alignment_length_aa: ", unique(width(aln_set))),
    "model: JTT + Inv + Gamma",
    paste0("bootstrap_replicates: ", opts$bootstrap),
    paste0("seed: ", opts$seed)
  ),
  con = summary_path
)

cat("Saved outputs in:", normalizePath(opts$out_dir), "\n")
cat(" -", alignment_path, "\n")
cat(" -", ml_tree_path, "\n")
cat(" -", ml_bs_tree_path, "\n")
cat(" -", pdf_path, "\n")
cat(" -", summary_path, "\n")
