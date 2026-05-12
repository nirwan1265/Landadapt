suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

base_gwas <- "/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt"
out_dir <- "Figs/Supplementary"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_order <- c("Arabidopsis", "Barley", "Rice", "Maize", "Sorghum")
species_to_crop <- c(
  Arabidopsis = "arabidopsis",
  Barley = "barley",
  Rice = "rice",
  Maize = "maize",
  Sorghum = "sorghum"
)

trait_cfg <- list(
  list(id = "pc1_worldclim", token = "PC1", label = "WorldClim PC1"),
  list(id = "pc2_worldclim", token = "PC2", label = "WorldClim PC2"),
  list(id = "pc3_worldclim", token = "PC3", label = "WorldClim PC3"),
  list(id = "ph", token = "ph", label = "Soil pH"),
  list(id = "soilN", token = "soilN", label = "Soil Nitrogen"),
  list(id = "am_rel_abundance_colonization", token = "am_rel_abundance_colonization", label = "AM Fungal Relative Abundance"),
  list(id = "am_roots_colonized", token = "am_roots_colonized", label = "AM Fungal Roots Colonized"),
  list(id = "aridity_index", token = "aridity_index", label = "Aridity Index")
)

make_path <- function(sp, token) {
  sp_lc <- species_to_crop[[sp]]
  prefix <- if (sp_lc == "arabidopsis") "AT" else sp_lc
  file.path(base_gwas, sp_lc, paste0(prefix, "_", token, ".txt"))
}

pick_col <- function(nms, candidates) {
  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NA_character_)
  nms[idx[[1]]]
}

build_qq_df <- function(pvals, species, max_points = 150000L) {
  p <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  n <- length(p)
  if (n < 10) return(NULL)
  p <- sort(p)
  idx <- seq_len(n)
  if (n > max_points) {
    keep <- unique(round(seq(1, n, length.out = max_points)))
    p <- p[keep]
    idx <- keep
  }
  expected <- -log10((idx - 0.5) / n)
  observed <- -log10(p)
  data.frame(species = species, expected = expected, observed = observed)
}

for (tr in trait_cfg) {
  qq_list <- list()
  for (sp in species_order) {
    f <- make_path(sp, tr$token)
    if (!file.exists(f)) {
      message("[WARN] Missing GWAS file: ", f)
      next
    }
    d <- readr::read_tsv(f, show_col_types = FALSE, progress = FALSE)
    p_col <- pick_col(names(d), c("p_wald", "P", "pvalue", "P.value", "p_lrt", "p_score"))
    if (is.na(p_col)) {
      message("[WARN] No p-value column in: ", f)
      next
    }
    pvals <- suppressWarnings(as.numeric(d[[p_col]]))
    qq <- build_qq_df(pvals, sp)
    if (!is.null(qq)) qq_list[[sp]] <- qq
  }

  if (length(qq_list) == 0) next
  qq_df <- bind_rows(qq_list) %>% mutate(species = factor(species, levels = species_order))

  max_xy <- max(c(qq_df$expected, qq_df$observed), na.rm = TRUE)
  p <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 0.2, alpha = 0.5, color = "black") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.35) +
    facet_wrap(~ species, ncol = 1, scales = "fixed") +
    coord_cartesian(xlim = c(0, max_xy), ylim = c(0, max_xy)) +
    labs(
      title = paste0("QQ Plots (", tr$label, " GWAS)"),
      x = expression(Expected~~-log[10](P)),
      y = expression(Observed~~-log[10](P))
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 12)
    )

  out <- file.path(out_dir, paste0("QQ_", tr$id, ".png"))
  ggsave(out, p, width = 8, height = 14, dpi = 300, bg = "white", device = ragg::agg_png)
  message("Saved: ", out)
}
