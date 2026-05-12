suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(tidyr)
})

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------
go_dir <- "results/GO"
out_dir <- file.path(go_dir, "common_pathways")
fig_dir <- "Figs/Fig3"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fractions <- c(0.01, 0.05) # top 1% and top 5%
top_n <- 50L
species_order <- c("arabidopsis", "barley", "rice", "maize", "sorghum")
drop_root_go <- c("GO:0008150", "GO:0003674", "GO:0005575")
species_colors <- c(
  arabidopsis = "#009E73",
  barley = "#0072B2",
  rice = "#CC79A7",
  maize = "#E69F00",
  sorghum = "#D55E00"
)

# ------------------------------------------------------------
# Load GO files
# ------------------------------------------------------------
go_files <- list.files(go_dir, pattern = "_GO_(BP|MF)\\.txt$", full.names = TRUE)
if (length(go_files) == 0) stop("No GO files found in ", go_dir)

read_go <- function(path) {
  nm <- basename(path)
  species <- tolower(str_extract(nm, "^[^_]+"))
  domain <- str_extract(nm, "(?<=_GO_)[^.]+")
  df <- read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  term_col <- names(df)[1]
  p_col <- names(df)[grepl("raw P-value", names(df), fixed = TRUE)][1]

  df %>%
    transmute(
      species = species,
      domain = domain,
      term_full = .data[[term_col]],
      p_raw = as.numeric(.data[[p_col]])
    ) %>%
    filter(!is.na(p_raw)) %>%
    mutate(
      go_id = str_extract(term_full, "GO:[0-9]{7}"),
      term = str_trim(str_remove(term_full, "\\s*\\(GO:[0-9]{7}\\)$"))
    ) %>%
    filter(!is.na(go_id), go_id != "")
}

go_all <- bind_rows(lapply(go_files, read_go))

select_top_n <- function(df, n_keep = 25L) {
  df %>% arrange(p_raw) %>% slice_head(n = n_keep)
}

make_outputs_for_top_n <- function(n_keep) {
  label <- paste0("top", n_keep, "terms")

  top_df <- go_all %>%
    group_by(species, domain) %>%
    group_modify(~ select_top_n(.x, n_keep)) %>%
    ungroup()

  top_df_noroot <- top_df %>% filter(!go_id %in% drop_root_go)

  common_df <- top_df_noroot %>%
    group_by(domain, go_id, term) %>%
    summarise(
      n_species = n_distinct(species),
      species_list = paste(sort(unique(species)), collapse = ";"),
      .groups = "drop"
    ) %>%
    filter(n_species >= 2) %>%
    arrange(domain, desc(n_species), go_id)

  common_long <- top_df_noroot %>%
    inner_join(common_df %>% select(domain, go_id, term, n_species), by = c("domain", "go_id", "term")) %>%
    mutate(logp = -log10(p_raw))

  common_wide <- common_long %>%
    select(domain, go_id, term, species, p_raw) %>%
    mutate(species = factor(species, levels = species_order)) %>%
    arrange(domain, go_id, species) %>%
    pivot_wider(names_from = species, values_from = p_raw)

  write.csv(top_df, file.path(out_dir, paste0("go_", label, "_selected.csv")), row.names = FALSE)
  write.csv(common_df, file.path(out_dir, paste0("go_", label, "_common_terms_ge2species.csv")), row.names = FALSE)
  write.csv(common_wide, file.path(out_dir, paste0("go_", label, "_common_terms_wide.csv")), row.names = FALSE)

  if (nrow(common_long) == 0) {
    p <- ggplot() +
      annotate("text", x = 1, y = 1, label = paste0("No non-root common GO terms for top ", n_keep, " terms"), size = 5) +
      xlim(0, 2) + ylim(0, 2) +
      theme_void()
  } else {
    term_levels <- common_df %>%
      mutate(term_lab = paste0(term, " (", go_id, ")")) %>%
      arrange(desc(n_species), term_lab) %>%
      pull(term_lab)

    plot_df <- common_long %>%
      mutate(
        species = factor(species, levels = species_order),
        term_lab = paste0(term, " (", go_id, ")"),
        term_lab = factor(term_lab, levels = rev(unique(term_levels)))
      )

    p <- ggplot(plot_df, aes(x = logp, y = term_lab, color = species)) +
      geom_point(size = 3.2, alpha = 0.95) +
      facet_wrap(~ domain, scales = "free_y") +
      scale_color_manual(values = species_colors, drop = FALSE) +
      labs(
        title = paste0("Common GO Pathways (Top ", n_keep, " terms by raw p-value per species)"),
        subtitle = "Only GO terms present in >=2 species (root terms removed)",
        x = expression(-log[10]("raw p-value")),
        y = "GO pathway / function",
        color = "Species"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
      )
  }

  fig_file <- file.path(fig_dir, paste0("GO_common_", label, "_pathways_dotplot.png"))
  ggsave(fig_file, p, width = 12, height = 6, dpi = 300, bg = "white", device = ragg::agg_png)

  message("Saved: ", fig_file)
  message("Common terms (>=2 species) using top ", n_keep, " terms: ", nrow(common_df))
}

make_outputs_for_top_n(top_n)
message("Done. Outputs in: ", out_dir)
