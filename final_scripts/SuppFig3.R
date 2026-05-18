suppressPackageStartupMessages({
  library(ggplot2)
})

base_dir <- "/Users/nirwantandukar/Documents/Github/Landadapt"
gwas_dir <- "/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt"
out_dir <- file.path(base_dir, "Figs", "Supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x = element_text(
      size = 24,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 24,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 24,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 24,
      color = "black"
    ),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 16
    ),
    plot.margin = margin(15, 15, 15, 15)
  )

species_info <- data.frame(
  species = c("Arabidopsis", "Barley", "Rice", "Maize", "Sorghum"),
  folder = c("arabidopsis", "barley", "rice", "maize", "sorghum"),
  prefix = c("AT", "barley", "rice", "maize", "sorghum"),
  color = c("#009E73", "#0072B2", "#CC79A7", "#E69F00", "#D55E00"),
  stringsAsFactors = FALSE
)

trait_info <- data.frame(
  trait = c("PC1", "PC2", "PC3", "pH", "soilN", "AM_rel", "AM_roots", "aridity"),
  suffix = c("PC1", "PC2", "PC3", "ph", "soilN", "am_rel_abundance_colonization", "am_roots_colonized", "aridity_index"),
  label = c(
    "WorldClim PC1",
    "WorldClim PC2",
    "WorldClim PC3",
    "Soil pH",
    "Soil Nitrogen",
    "AM Fungal Relative Abundance",
    "AM Fungal Roots Colonized",
    "Aridity Index"
  ),
  stringsAsFactors = FALSE
)

read_gwas <- function(species_row, trait_row) {
  filename <- paste0(species_row[["prefix"]], "_", trait_row[["suffix"]], ".txt")
  path <- file.path(gwas_dir, species_row[["folder"]], filename)
  if (!file.exists(path)) {
    stop("Missing GWAS file: ", path)
  }

  header_fields <- strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
  p_candidates <- c("p_wald", "P.value", "pvalue", "p_value", "P_VALUE")
  p_idx <- match(TRUE, header_fields %in% p_candidates)
  if (is.na(p_idx)) {
    stop("No supported p-value column found in: ", path)
  }
  col_classes <- rep("NULL", length(header_fields))
  col_classes[p_idx] <- "numeric"
  dat <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    colClasses = col_classes
  )

  p_col <- names(dat)[1]
  pvals <- suppressWarnings(as.numeric(dat[[p_col]]))
  pvals <- pvals[is.finite(pvals) & !is.na(pvals) & pvals > 0 & pvals <= 1]
  pvals <- sort(pvals, decreasing = FALSE)

  n <- length(pvals)
  if (n == 0) {
    stop("No valid p-values found in: ", path)
  }

  max_points <- 5000L
  if (n > max_points) {
    keep_idx <- unique(round(seq(1, n, length.out = max_points)))
    pvals <- pvals[keep_idx]
    n <- length(pvals)
  }

  data.frame(
    species = species_row[["species"]],
    trait = trait_row[["trait"]],
    trait_label = trait_row[["label"]],
    expected = -log10(ppoints(n)),
    observed = -log10(pvals),
    stringsAsFactors = FALSE
  )
}

qq_df <- do.call(
  rbind,
  lapply(seq_len(nrow(trait_info)), function(i) {
    do.call(
      rbind,
      lapply(seq_len(nrow(species_info)), function(j) {
        read_gwas(species_info[j, ], trait_info[i, ])
      })
    )
  })
)

qq_df$trait_label <- factor(qq_df$trait_label, levels = trait_info$label)
qq_df$species <- factor(qq_df$species, levels = species_info$species)

x_max <- ceiling(max(qq_df$expected))
y_max <- ceiling(max(qq_df$observed))

p <- ggplot(qq_df, aes(x = expected, y = observed, color = species)) +
  geom_segment(
    x = 0, y = 0, xend = x_max, yend = x_max,
    linetype = "dashed", linewidth = 0.5, color = "grey45"
  ) +
  geom_line(linewidth = 0.8, alpha = 0.95) +
  facet_wrap(~trait_label, ncol = 2, scales = "fixed") +
  scale_color_manual(values = setNames(species_info$color, species_info$species)) +
  scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "SuppFig3: QQ plots across traits (species overlay)",
    x = expression("Expected " * -log[10] * "(P)"),
    y = expression("Observed " * -log[10] * "(P)")
  ) +
  plot_theme +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

png_path <- file.path(out_dir, "SuppFig3_QQ_all_traits_species_overlay.png")
pdf_path <- file.path(out_dir, "SuppFig3_QQ_all_traits_species_overlay.pdf")

ggsave(png_path, p, width = 14, height = 18, dpi = 300, bg = "white")
ggsave(pdf_path, p, width = 14, height = 18, bg = "white")

cat("WROTE", png_path, "\n")
cat("WROTE", pdf_path, "\n")
