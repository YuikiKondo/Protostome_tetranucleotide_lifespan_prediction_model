setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/04_mechanical_index_and_lifespan_correlation")

# =======================
# Setup
# =======================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(tidyr)
})

# -----------------------
# Inputs
# -----------------------
input_file <- "filtered_species_with_all_indexes_di_tri_tetra_combined_Github_check.csv"
out_dir    <- "scatterplots_ln_lifespan_by_region_batch_with_stats"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

indices <- c(
  "DNA_bendability","X_Displacement","Y_Displacement","Inclination","Tip",
  "Shear","Stretch","Stagger","Buckle","Propel","Opening","Alpha","Beta","Gamma",
  "Delta","Epsilon","Zeta","Chi","Phase","Amplitude","Hydrogen_bond",
  "Stacking_energy","Solvation","Shift","Slide","Rise","Tilt","Roll","Twist",
  "A_philicity_energy"
)

# -----------------------
# Load & prepare data
# -----------------------
df <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    ln_Average_lifespan_days = ifelse(Average_lifespan_days > 0,
                                      log(Average_lifespan_days), NA_real_),
    Region = case_when(
      Region == "upstream2"   ~ "Upstream",
      Region == "downstream2" ~ "Downstream",
      Region == "genome"      ~ "Genome",
      Region == "exon"        ~ "Exon",
      Region == "intron"      ~ "Intron",
      Region %in% c("upstream1", "downstream1") ~ NA_character_,
      TRUE ~ str_to_title(Region)
    )
  ) %>%
  filter(!is.na(Region), !is.na(ln_Average_lifespan_days))

region_levels <- c("Genome", "Upstream", "Exon", "Intron", "Downstream")
df$Region <- factor(df$Region, levels = region_levels)

# -----------------------
# Compute Pearson r and global FDR across ALL panels
# -----------------------
# long format across all 30 indices
df_long_all <- df %>%
  select(ln_Average_lifespan_days, Region, all_of(indices)) %>%
  pivot_longer(cols = all_of(indices),
               names_to = "Index",
               values_to = "Value")

# correlation per (Index, Region)
stats_all <- df_long_all %>%
  group_by(Index, Region) %>%
  reframe({
    d <- drop_na(cur_data(), ln_Average_lifespan_days, Value)
    n <- nrow(d)
    if (n >= 3 && sd(d$ln_Average_lifespan_days) > 0 && sd(d$Value) > 0) {
      ct <- suppressWarnings(cor.test(d$ln_Average_lifespan_days, d$Value, method = "pearson"))
      tibble(
        n = n,
        r = unname(ct$estimate),
        p = ct$p.value
      )
    } else {
      tibble(n = n, r = NA_real_, p = NA_real_)
    }
  }) %>%
  ungroup()

# Global BH/FDR across all Index×Region tests
stats_all <- stats_all %>%
  mutate(q = p.adjust(p, method = "BH"),
         Index_pretty = str_replace_all(Index, "_", " "))

# helper: label text formatting
fmt_r <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))
fmt_q <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
stats_all <- stats_all %>%
  mutate(label = paste0("r = ", fmt_r(r), "\nq = ", fmt_q(q)))

# -----------------------
# Batch plotting: 6 indices per PNG, with stats
# -----------------------
batches <- split(indices, ceiling(seq_along(indices)/6))

walk2(batches, seq_along(batches), function(batch_indices, batch_num) {

  dat_long <- df %>%
    select(ln_Average_lifespan_days, Region, all_of(batch_indices)) %>%
    pivot_longer(cols = all_of(batch_indices),
                 names_to = "Index",
                 values_to = "Value") %>%
    mutate(Index_pretty = str_replace_all(Index, "_", " "))

  # attach stats (precomputed globally) for these indices
  stats_batch <- stats_all %>%
    filter(Index %in% batch_indices) %>%
    select(Index, Index_pretty, Region, label)

  p <- ggplot(dat_long, aes(x = ln_Average_lifespan_days, y = Value)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    facet_grid(rows = vars(Index_pretty), cols = vars(Region),
               scales = "free_y", switch = "y") +
    labs(x = "ln(lifespan days)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text.x = element_text(face = "bold"),
      strip.text.y.left = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.placement = "outside",
      panel.spacing = unit(0.8, "lines")
    ) +
    # Add per-panel labels at top-left (works with free y using -Inf/Inf)
    # Add per-panel labels at top-left (works with free y using -Inf/Inf)
    geom_text(
      data = stats_batch %>%
        mutate(
          label_color = ifelse(as.numeric(str_extract(label, "(?<=q = )[0-9\\.eE\\-]+")) < 0.05, "red", "black"),
          label_font  = ifelse(as.numeric(str_extract(label, "(?<=q = )[0-9\\.eE\\-]+")) < 0.05, "bold", "plain")
        ),
      mapping = aes(x = -Inf, y = Inf, label = label, color = label_color, fontface = label_font),
      hjust = -0.1, vjust = 1.1, size = 3, show.legend = FALSE
    ) +
    scale_color_identity()


  ggsave(file.path(out_dir, paste0("scatter_batch_", batch_num, ".png")),
         p, width = 14, height = 12, dpi = 300)
})

message("Done. PNGs (with r & FDR q on each panel) are in: ", normalizePath(out_dir))



# -----------------------
# Save stats table as CSV
# -----------------------
stats_outfile <- file.path(out_dir, "pearson_correlation_results_DNA_mechanics_and_lifespan_global_FDR.csv")

stats_all %>%
  select(Index, Region, n, r, p, FDR_adjusted_p = q) %>%
  arrange(Index, Region) %>%
  write_csv(stats_outfile)

message("Correlation results saved to: ", normalizePath(stats_outfile))
