suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(forcats)
  library(patchwork)
})

# -----------------------------
# Config
# -----------------------------
go_dir <- "results/GO"
out_dir <- file.path(go_dir, "top1pct")
fig_dir <- "Figs/Fig3"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

top_fraction <- 0.01
species_order <- c("arabidopsis", "barley", "rice", "maize", "sorghum")

# -----------------------------
# Helpers
# -----------------------------
parse_go_file <- function(path) {
  nm <- basename(path)
  sp <- str_extract(nm, "^[^_]+")
  dom <- str_extract(nm, "(?<=_GO_)[^.]+")

  df <- read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  term_col <- names(df)[1]
  p_col <- names(df)[grepl("raw P-value", names(df), fixed = TRUE)]
  if (length(p_col) == 0) stop("No raw P-value column in: ", path)
  p_col <- p_col[1]

  out <- df %>%
    transmute(
      species = tolower(sp),
      domain = dom,
      term_full = .data[[term_col]],
      p_raw = as.numeric(.data[[p_col]])
    ) %>%
    filter(!is.na(p_raw)) %>%
    mutate(
      go_id = str_extract(term_full, "GO:[0-9]{7}"),
      term = str_trim(str_remove(term_full, "\\s*\\(GO:[0-9]{7}\\)$"))
    ) %>%
    filter(!is.na(go_id), go_id != "")

  out
}

select_top_fraction <- function(df, frac = 0.01) {
  n_total <- nrow(df)
  n_keep <- max(1L, floor(n_total * frac))
  out <- df %>% arrange(p_raw) %>% slice_head(n = n_keep)
  list(data = out, n_total = n_total, n_keep = n_keep)
}

# -----------------------------
# Load + select top 1%
# -----------------------------
go_files <- list.files(go_dir, pattern = "_GO_(BP|MF)\\.txt$", full.names = TRUE)
if (length(go_files) == 0) stop("No GO files found in: ", go_dir)

parsed <- lapply(go_files, parse_go_file)
top_list <- lapply(parsed, select_top_fraction, frac = top_fraction)
top_df <- bind_rows(lapply(top_list, `[[`, "data"))

summary_df <- bind_rows(lapply(top_list, function(x) {
  d <- x$data
  data.frame(
    species = unique(d$species),
    domain = unique(d$domain),
    n_total_terms = x$n_total,
    n_top_selected = x$n_keep,
    stringsAsFactors = FALSE
  )
}))

write.csv(top_df, file.path(out_dir, "go_top1pct_selected_all.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "go_top1pct_selection_summary.csv"), row.names = FALSE)

# Save per species-domain tables
split_keys <- unique(top_df %>% select(species, domain))
for (i in seq_len(nrow(split_keys))) {
  sp <- split_keys$species[i]
  dm <- split_keys$domain[i]
  out <- top_df %>% filter(species == sp, domain == dm) %>% arrange(p_raw)
  write.csv(out, file.path(out_dir, paste0(sp, "_", dm, "_top1pct.csv")), row.names = FALSE)
}

# -----------------------------
# Common GO terms (>=2 species) within each domain
# -----------------------------
common_all <- top_df %>%
  group_by(domain, go_id, term) %>%
  summarise(
    n_species = n_distinct(species),
    species_list = paste(sort(unique(species)), collapse = ";"),
    .groups = "drop"
  ) %>%
  filter(n_species >= 2) %>%
  arrange(domain, desc(n_species), go_id)

write.csv(common_all, file.path(out_dir, "go_top1pct_common_terms_ge2species.csv"), row.names = FALSE)
write.csv(common_all %>% filter(domain == "BP"), file.path(out_dir, "go_top1pct_common_BP_ge2species.csv"), row.names = FALSE)
write.csv(common_all %>% filter(domain == "MF"), file.path(out_dir, "go_top1pct_common_MF_ge2species.csv"), row.names = FALSE)

# -----------------------------
# Figure 1: common GO heatmaps (BP/MF)
# -----------------------------
plot_domain <- function(domain_name) {
  d_common <- common_all %>% filter(domain == domain_name)

  if (nrow(d_common) == 0) {
    return(
      ggplot() +
        annotate("text", x = 1, y = 1, label = paste0("No common ", domain_name, "\nGO terms at top 1%"), size = 5) +
        xlim(0, 2) + ylim(0, 2) +
        theme_void() +
        ggtitle(paste0(domain_name, " (Top 1% Raw p-value)"))
    )
  }

  long_df <- top_df %>%
    filter(domain == domain_name) %>%
    select(species, domain, go_id, term, p_raw) %>%
    inner_join(d_common %>% select(domain, go_id), by = c("domain", "go_id")) %>%
    mutate(logp = -log10(p_raw))

  species_levels <- species_order[species_order %in% unique(long_df$species)]
  term_levels <- d_common %>%
    arrange(desc(n_species), go_id) %>%
    pull(term) %>%
    unique()

  plot_df <- expand_grid(
    species = species_levels,
    term = term_levels
  ) %>%
    left_join(
      long_df %>% select(species, term, logp),
      by = c("species", "term")
    )

  ggplot(plot_df, aes(x = species, y = term, fill = logp)) +
    geom_tile(color = "grey90") +
    scale_fill_gradient(low = "white", high = "black", na.value = "grey95", name = expression(-log[10](p))) +
    labs(
      title = paste0(domain_name, " Common GO Terms (Top 1%)"),
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

p_bp <- plot_domain("BP")
p_mf <- plot_domain("MF")
fig <- p_bp / p_mf + plot_annotation(tag_levels = "A")

fig_out <- file.path(fig_dir, "GO_common_top1pct.png")
ggsave(fig_out, fig, width = 12, height = 10, dpi = 300, bg = "white", device = ragg::agg_png)

# -----------------------------
# Figure 2: common GO dotplot (x = -log10 raw p, y = function)
# -----------------------------
common_long <- top_df %>%
  inner_join(common_all %>% select(domain, go_id, term), by = c("domain", "go_id", "term")) %>%
  group_by(species, domain, go_id, term) %>%
  summarise(p_raw = min(p_raw, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    logp = -log10(p_raw),
    species = factor(species, levels = species_order),
    term_y = paste0(term, " (", go_id, ")")
  )

if (nrow(common_long) > 0) {
  p_dot <- ggplot(common_long, aes(x = logp, y = fct_reorder(term_y, logp, .fun = max), color = species)) +
    geom_point(size = 3, alpha = 0.95) +
    facet_wrap(~ domain, scales = "free_y") +
    labs(
      title = "Common GO Terms Across Species (Top 1%)",
      subtitle = "Terms present in >=2 species",
      x = expression(-log[10]("raw p-value")),
      y = "GO function",
      color = "Species"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
} else {
  p_dot <- ggplot() +
    annotate("text", x = 1, y = 1, label = "No common GO terms at top 1%", size = 5) +
    xlim(0, 2) + ylim(0, 2) +
    theme_void()
}

fig_dot_out <- file.path(fig_dir, "GO_common_top1pct_dotplot.png")
ggsave(fig_dot_out, p_dot, width = 11, height = 5, dpi = 300, bg = "white", device = ragg::agg_png)

message("Saved GO top1% tables to: ", out_dir)
message("Saved figure: ", fig_out)
message("Saved figure: ", fig_dot_out)
