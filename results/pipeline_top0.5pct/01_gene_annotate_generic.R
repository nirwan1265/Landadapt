#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 01_gene_annotate_generic.R <gwas_file> <gff3_file> <output_file> [window_bp]", call. = FALSE)
}

gwas_file <- args[[1]]
gff3_file <- args[[2]]
out_file <- args[[3]]
window_bp <- if (length(args) >= 4) as.integer(args[[4]]) else 25000L

if (!file.exists(gwas_file)) stop("GWAS file not found: ", gwas_file, call. = FALSE)
if (!file.exists(gff3_file)) stop("GFF3 file not found: ", gff3_file, call. = FALSE)
if (is.na(window_bp) || window_bp < 0) stop("window_bp must be >=0", call. = FALSE)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

normalize_chr <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- sub("^9311_", "", x)
  x <- sub("^chromosome", "", x)
  x <- sub("^chrom", "", x)
  x <- sub("^chr", "", x)
  # trim leading zeros only for pure numeric IDs
  is_num <- grepl("^[0-9]+$", x)
  x[is_num] <- as.character(as.integer(x[is_num]))
  x
}

extract_attr <- function(attr, key) {
  pat <- paste0(".*(?:^|;)", key, "=([^;]+).*")
  out <- sub(pat, "\\1", attr, perl = TRUE)
  out[out == attr] <- NA_character_
  out
}

message("[1/4] Read GFF3 genes: ", gff3_file)
gff <- fread(
  cmd = paste("grep -v '^#'", shQuote(gff3_file)),
  sep = "\t", header = FALSE, quote = "", fill = TRUE, showProgress = FALSE
)
if (ncol(gff) < 9) stop("Malformed GFF3 (need >=9 cols): ", gff3_file, call. = FALSE)
setnames(gff, paste0("V", seq_len(ncol(gff))))

genes <- gff[V3 == "gene", .(
  chr_raw = V1,
  gene_start = as.integer(V4),
  gene_end = as.integer(V5),
  attrs = V9
)]
if (nrow(genes) == 0) stop("No gene features in GFF3: ", gff3_file, call. = FALSE)

genes[, chr_norm := normalize_chr(chr_raw)]
# prefer Name, then gene_id, then ID stripped
name <- extract_attr(genes$attrs, "Name")
gene_id_attr <- extract_attr(genes$attrs, "gene_id")
id <- extract_attr(genes$attrs, "ID")
id <- sub("^gene:", "", id)
final_gene <- fifelse(!is.na(name) & name != "", name,
                 fifelse(!is.na(gene_id_attr) & gene_id_attr != "", gene_id_attr, id))

genes[, gene_name := utils::URLdecode(final_gene)]
genes[, gene_id := utils::URLdecode(id)]
genes[, window_start := pmax(1L, gene_start - window_bp)]
genes[, window_end := gene_end + window_bp]

gene_windows <- genes[, .(chr_norm, window_start, window_end, gene_name, gene_id, gene_start, gene_end)]
setkey(gene_windows, chr_norm, window_start, window_end)

message("[2/4] Read GWAS: ", gwas_file)
gwas <- fread(gwas_file, showProgress = FALSE)

pick_col <- function(cands, nms) {
  hit <- cands[cands %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

chr_col <- pick_col(c("chr", "CHR", "chrom", "Chrom", "CHROM", "chromosome", "Chromosome", "Chr"), names(gwas))
pos_col <- pick_col(c("ps", "BP", "bp", "pos", "POS", "position", "Position", "Pos"), names(gwas))
if (is.na(chr_col) || is.na(pos_col)) {
  stop("Could not detect chr/pos columns in GWAS: ", gwas_file,
       "\ncolumns=", paste(names(gwas), collapse = ", "), call. = FALSE)
}

gwas[, row_id := .I]
gwas[, chr_norm := normalize_chr(get(chr_col))]
gwas[, pos_i := as.integer(get(pos_col))]

message("[3/4] Overlap SNPs with gene windows (+/-", window_bp, " bp)")
chr_vals <- unique(gwas$chr_norm)
res <- vector("list", length(chr_vals))
ri <- 1L
for (cv in chr_vals) {
  snp_chr <- gwas[chr_norm == cv & !is.na(pos_i), .(row_id, pos_i, start = pos_i, end = pos_i)]
  if (nrow(snp_chr) == 0) next
  genes_chr <- gene_windows[chr_norm == cv]
  if (nrow(genes_chr) == 0) next

  snp_chr[, chr_norm := cv]
  setkey(snp_chr, chr_norm, start, end)

  ov <- foverlaps(
    x = snp_chr,
    y = genes_chr,
    by.x = c("chr_norm", "start", "end"),
    by.y = c("chr_norm", "window_start", "window_end"),
    nomatch = 0L
  )
  if (nrow(ov) == 0) next

  ov[, dist_bp := fifelse(pos_i < gene_start, gene_start - pos_i,
                   fifelse(pos_i > gene_end, pos_i - gene_end, 0L))]

  ann <- ov[order(row_id, dist_bp, gene_name), .(
    n_genes_25kb = uniqueN(gene_name),
    genes_in_25kb = paste(unique(gene_name), collapse = ";"),
    closest_gene = gene_name[1L],
    closest_gene_id = gene_id[1L],
    closest_gene_distance_bp = dist_bp[1L]
  ), by = row_id]

  res[[ri]] <- ann
  ri <- ri + 1L
}

ann_all <- rbindlist(res, use.names = TRUE, fill = TRUE)
if (nrow(ann_all) == 0) {
  ann_all <- data.table(
    row_id = integer(),
    n_genes_25kb = integer(),
    genes_in_25kb = character(),
    closest_gene = character(),
    closest_gene_id = character(),
    closest_gene_distance_bp = integer()
  )
}

message("[4/4] Write annotated file: ", out_file)
out <- merge(gwas, ann_all, by = "row_id", all.x = TRUE, sort = FALSE)
out[is.na(n_genes_25kb), n_genes_25kb := 0L]
# drop helper cols
out[, c("row_id", "chr_norm", "pos_i") := NULL]
fwrite(out, out_file, sep = "\t", quote = FALSE, na = "NA")

message("DONE")
