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

species_colors <- c(
  Arabidopsis = "#009E73",
  Barley = "#0072B2",
  Rice = "#CC79A7",
  Maize = "#E69F00",
  Sorghum = "#D55E00"
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

build_qq_df <- function(pvals, species, trait_label, max_points = 120000L) {
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
  data.frame(
    species = species,
    trait = trait_label,
    expected = expected,
    observed = observed
  )
}

qq_all <- list()
for (tr in trait_cfg) {
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
    qq <- build_qq_df(pvals, sp, tr$label)
    if (!is.null(qq)) qq_all[[paste(tr$id, sp, sep = "_")]] <- qq
  }
}

qq_df <- bind_rows(qq_all)
if (nrow(qq_df) == 0) stop("No QQ data available for combined figure.")

qq_df <- qq_df %>%
  mutate(
    species = factor(species, levels = species_order),
    trait = factor(trait, levels = vapply(trait_cfg, `[[`, character(1), "label"))
  )

diag_df <- qq_df %>%
  group_by(trait) %>%
  summarise(max_xy = max(c(expected, observed), na.rm = TRUE), .groups = "drop")

p <- ggplot(qq_df, aes(x = expected, y = observed, color = species, group = species)) +
  geom_line(linewidth = 0.35, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, color = "grey20", linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  scale_color_manual(values = species_colors, drop = FALSE) +
  labs(
    title = "SuppFig3: QQ Plots Across Traits (Species Overlay)",
    x = expression(Expected~~-log[10](P)),
    y = expression(Observed~~-log[10](P)),
    color = "Species"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12)
  )

out_png <- file.path(out_dir, "SuppFig3_QQ_all_traits_species_overlay.png")
out_pdf <- file.path(out_dir, "SuppFig3_QQ_all_traits_species_overlay.pdf")

ggsave(out_png, p, width = 14, height = 18, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(out_pdf, p, width = 14, height = 18, bg = "white", device = cairo_pdf)

message("Saved: ", out_png)
message("Saved: ", out_pdf)
