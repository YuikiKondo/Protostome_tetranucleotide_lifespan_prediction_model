# =======================
# Setup
# =======================
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/05_mechanical_index_parameter_and_weight_correlation")

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tools)

# Input (by-region) produced by your Python step
infile <- "tetranucleotide_mechanics_full_zscore_with_weights_by_region.csv"
if (!file.exists(infile)) {
  stop("Input file not found: tetranucleotide_mechanics_full_zscore_with_weights_by_region.csv")
}
df <- read_csv(infile, show_col_types = FALSE)
stopifnot("Motif" %in% names(df))

# =======================
# Identify columns
# =======================
# All weight columns (per-region)
weight_cols_all <- grep("^Average_Weight_", names(df), value = TRUE)
# Exclude the across-region mean if present
region_cols <- setdiff(weight_cols_all, "Average_Weight_Tetranuc")
if (length(region_cols) == 0) stop("No region-specific weight columns found.")

# Mechanics columns: *_z, numeric, and not weight columns
exclude_cols <- c("Motif", weight_cols_all)
mech_cols <- setdiff(names(df), exclude_cols)
mech_cols <- mech_cols[grepl("_z$", mech_cols)]
mech_cols <- mech_cols[sapply(df[mech_cols], is.numeric)]
if (length(mech_cols) == 0) stop("No mechanics *_z numeric columns found.")

# =======================
# Helpers
# =======================
spearman_test <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  x2 <- x[ok]; y2 <- y[ok]
  if (length(x2) < 3L) return(list(rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x2, y2, method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value)
}

fmt_num <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
fmt_p   <- function(p) ifelse(is.na(p), "NA", format.pval(p, digits = 2, eps = 1e-3))

pretty_mech_label <- function(colname) {
  base <- sub("_z$", "", colname)
  pretty <- gsub("_", " ", base)
  paste0(pretty, " (Z-scored)")
}

region_key <- function(col) sub("^Average_Weight_", "", col)
pretty_region <- function(reg) toTitleCase(reg)

# enforce your preferred region order
desired_regions <- c("exon", "intron", "upstream", "downstream", "genome")

# =======================
# 1) Compute Spearman across ALL (index × region) tests
# =======================
results_list <- list()
for (rc in region_cols) {
  reg <- region_key(rc)
  for (idx in mech_cols) {
    st <- spearman_test(df[[rc]], df[[idx]])
    results_list[[length(results_list) + 1L]] <- tibble::tibble(
      Region        = reg,
      Weight_Column = rc,
      Index         = idx,
      Spearman_rho  = st$rho,
      p_value       = st$p
    )
  }
}
all_results <- dplyr::bind_rows(results_list)

# =======================
# 2) Global FDR (BH) across ALL tests
# =======================
all_results <- all_results %>%
  mutate(p_value_FDR_corrected_GLOBAL = p.adjust(p_value, method = "BH"))

# Save combined table
outroot <- "mechanics_scatterplots_fdr_by_region_GLOBAL"
dir.create(outroot, showWarnings = FALSE)
readr::write_csv(all_results,
                 file.path(outroot, "spearman_p_and_GLOBAL_FDR_summary_ALL_REGIONS.csv"))

# =======================
# 3) Long data for combined-panel plots (one row per mechanics index)
# =======================
long_df <- df %>%
  pivot_longer(
    cols = all_of(region_cols),
    names_to = "Weight_Column",
    values_to = "Average_Weight"
  ) %>%
  mutate(
    Region = region_key(Weight_Column),
    Region_pretty = factor(pretty_region(Region),
                           levels = pretty_region(desired_regions)[desired_regions %in% Region])
  )

# =======================
# 4) Plotting: for each mechanics index, make a single PNG with 5 facets (regions)
# =======================
combined_dir <- file.path(outroot, "combined_rows_per_index")
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

SIG_Q <- 0.05

for (idx in mech_cols) {
  # metrics for labels (this already carries global FDR)
  subres <- all_results %>%
    filter(Index == idx) %>%
    mutate(
      Region_pretty = factor(pretty_region(Region),
                             levels = pretty_region(desired_regions)[desired_regions %in% Region])
    )

  # label table (one row per facet)
  lab_tbl <- subres %>%
      transmute(
        Region_pretty = factor(pretty_region(Region),
          levels = pretty_region(desired_regions)[desired_regions %in% Region]
        ),
        label = paste0(
          "Spearman's rho = ", fmt_num(Spearman_rho),
          "\n", "p-value = ", fmt_p(p_value),
          "\n", "FDR p (global) = ", fmt_p(p_value_FDR_corrected_GLOBAL)
        ),
        is_sig = p_value_FDR_corrected_GLOBAL < SIG_Q,
        x = Inf, y = Inf, hjust = 1.05, vjust = 1.2
      )


  # plotting data for this index
  plot_df <- long_df %>%
    filter(Region %in% subres$Region)

  p <- ggplot(plot_df, aes(x = Average_Weight, y = .data[[idx]])) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "blue") +
  facet_wrap(~ Region_pretty, nrow = 1, scales = "free") +
  # annotation with padding
  # Significant: bold + red
  geom_text(
    data = dplyr::filter(lab_tbl, is_sig),
    mapping = aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE, size = 3, color = "red", fontface = "bold"
  ) +
  # Non-significant: black + plain
  geom_text(
    data = dplyr::filter(lab_tbl, !is_sig),
    mapping = aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE, size = 3, color = "black", fontface = "plain"
  ) +
  labs(
    x = "Average weight by region",
    y = pretty_mech_label(idx),
    title = pretty_mech_label(idx)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(t = 5, r = 20, b = 20, l = 20)  # more margin space
  )


  # width ~ 5 panels * 3.2in each; adjust as you like
  outfile <- file.path(
    combined_dir,
    paste0("scatter_ROW_", idx, "_vs_AverageWeight_GLOBALFDR.png")
  )
  ggsave(outfile, plot = p, width = 16, height = 4.8, dpi = 300)
}

message("✅ Done. Wrote combined (5-panel) PNGs per mechanics index to: ", combined_dir)




#### batched png #####

# =======================
# 4) Plotting in batches: 6 mechanics per PNG × 5 PNGs
# =======================
combined_dir <- file.path(outroot, "combined_rows_per_index_batched")
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

# Use the SAME pretty label function (with "(Z-scored)") everywhere
pretty_mech_strip <- pretty_mech_label  # ensure consistency

annot_all <- all_results %>%
  mutate(
    Region_pretty = factor(pretty_region(Region),
      levels = pretty_region(desired_regions)[desired_regions %in% Region]
    ),
    Index_strip = pretty_mech_strip(Index),
    label = paste0(
      "Spearman's rho = ", fmt_num(Spearman_rho),
      "\n", "p-value = ", fmt_p(p_value),
      "\n", "FDR p (global) = ", fmt_p(p_value_FDR_corrected_GLOBAL)
    ),
    is_sig = p_value_FDR_corrected_GLOBAL < SIG_Q,   # keep significance flag
    x = Inf, y = Inf, hjust = 1.02, vjust = 1.08
  )

chunk_size <- 6L
chunks <- split(mech_cols, ceiling(seq_along(mech_cols) / chunk_size))

for (i in seq_along(chunks)) {
  idx_set <- chunks[[i]]

  # 1) Data to plot
  plot_df <- long_df %>%
    select(Motif, Weight_Column, Average_Weight, Region, Region_pretty, all_of(idx_set)) %>%
    pivot_longer(
      cols = all_of(idx_set),
      names_to = "Index",
      values_to = "MechValue"
    ) %>%
    mutate(
      Index_strip = pretty_mech_strip(Index),   # full label with (Z-scored)
      Index_strip = factor(Index_strip, levels = pretty_mech_strip(idx_set)),
      Region_pretty = factor(Region_pretty,
        levels = pretty_region(desired_regions)[desired_regions %in% Region]
      )
    )

  # 2) Labels (must include is_sig)
  lab_tbl <- annot_all %>%
    filter(Index %in% idx_set) %>%
    mutate(Index_strip = factor(Index_strip, levels = pretty_mech_strip(idx_set))) %>%
    select(Region_pretty, Index_strip, label, is_sig, x, y, hjust, vjust)

  # 3) Plot
  p <- ggplot(plot_df, aes(x = Average_Weight, y = MechValue)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "blue") +
    facet_grid(
      rows = vars(Index_strip),
      cols = vars(Region_pretty),
      scales = "free_y",
      switch = "y"
    ) +
    # Significant labels: bold + red
    geom_text(
      data = dplyr::filter(lab_tbl, is_sig),
      aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE, size = 2.8, color = "red", fontface = "bold"
    ) +
    # Non-significant labels: black + plain
    geom_text(
      data = dplyr::filter(lab_tbl, !is_sig),
      aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE, size = 2.8, color = "black", fontface = "plain"
    ) +
    labs(x = "Average weight by region") +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 90, face = "bold", margin = margin(r = 6)),
      strip.text.x = element_text(face = "bold"),
      axis.title.y = element_blank(),
      plot.margin = margin(t = 8, r = 30, b = 12, l = 60)
    )

  outfile <- file.path(
    combined_dir,
    sprintf("scatter_BATCH_%02d_of_%02d_6indices_x_5regions_GLOBALFDR.png", i, length(chunks))
  )
  ggsave(outfile, plot = p, width = 16, height = 18, dpi = 300)
}


message("✅ Done. Wrote ", length(chunks),
        " batched PNGs (6 mechanics per file) to: ", combined_dir)

