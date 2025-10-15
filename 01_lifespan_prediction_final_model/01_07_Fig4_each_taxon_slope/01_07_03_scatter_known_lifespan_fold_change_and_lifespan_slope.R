# ====== Set working directory ======
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_final_model/01_07_Fig4_each_taxon_slope")

# ===================== Setup =====================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  # ggrepel is optional for labels; uncomment if you want labels later
  # library(ggrepel)
})


# ===================== 1) Locate files =====================
# Matches the CSVs you saved in your script
level_files <- list.files(pattern = "^slope_by_(phylum|class|order|family|genus)_ln_years\\.csv$", ignore.case = TRUE)
proto_file  <- "slope_all_protostomia_ln_years.csv"
files_to_read <- c(level_files, if (file.exists(proto_file)) proto_file else character())

if (length(files_to_read) == 0) {
  stop("No slope CSV files found in the working directory.")
}

# ===================== 2) Helper to parse Level from filename =====================
parse_level <- function(filename) {
  fn <- tolower(filename)
  if (str_detect(fn, "protostomia")) return("Protostomia")
  m <- str_match(fn, "^slope_by_(phylum|class|order|family|genus)_ln_years\\.csv$")
  if (!is.na(m[1,2])) {
    # Capitalize first letter
    return(str_to_title(m[1,2]))
  }
  NA_character_
}

# ===================== 3) Read, add Level, combine =====================
slope_df_list <- lapply(files_to_read, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  df$Level <- parse_level(f)
  df
})

slope_all <- bind_rows(slope_df_list)

# Sanity check required columns
required_cols <- c("slope", "Max_minus_Min_Known_Lifespan", "Level")
missing_cols  <- setdiff(required_cols, colnames(slope_all))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns in combined data: %s",
               paste(missing_cols, collapse = ", ")))
}

# Optional: ensure the taxon name column exists (for labels if wanted).
# Your per-level files have the focal grouping column named by the Level (e.g., 'Phylum', 'Class', ...).
# We'll create a unified 'Taxon' column by picking whichever of those exists on each row.
taxon_cols <- c("Phylum","Class","Order","Family","Genus")
slope_all <- slope_all %>%
  mutate(
    Taxon = coalesce(!!!rlang::syms(taxon_cols))  # first non-NA among those columns
  )

# Order legend nicely
slope_all <- slope_all %>%
  mutate(Level = factor(Level, levels = c("Protostomia","Phylum","Class","Order","Family","Genus")))

# ===================== 4) Plot: Max_minus_Min_Known_Lifespan vs slope =====================
# ===================== 4) Plot: Max_minus_Min_Known_Lifespan vs slope =====================
# --- 1) Compute median slope per Level ---
median_slope_by_level <- slope_all %>%
  dplyr::group_by(Level) %>%
  dplyr::summarise(median_slope = median(slope, na.rm = TRUE), .groups = "drop")

# --- 2) Run one-tailed Wilcoxon tests per Level (median slope > 0) ---
slope_wilcoxon_tests <- slope_all %>%
  dplyr::group_by(Level) %>%
  dplyr::group_modify(~{
    x <- .x$slope[is.finite(.x$slope)]
    n <- length(x)
    if (n > 2) {
      wt <- tryCatch(
        stats::wilcox.test(x, mu = 0, alternative = "greater"),
        error = function(e) NULL
      )
      tibble::tibble(
        n = n,
        wilcox_p = ifelse(is.null(wt), NA_real_, wt$p.value)
      )
    } else {
      tibble::tibble(n = n, wilcox_p = NA_real_)
    }
  }) %>%
  dplyr::ungroup()

# --- 3) Combine results into one legend label (median, n, p) ---
legend_info <- median_slope_by_level %>%
  dplyr::left_join(slope_wilcoxon_tests, by = "Level") %>%
  dplyr::mutate(
    label = sprintf(
      "%s (median slope = %.2f, n = %d, p = %.3g)",
      as.character(Level), median_slope, n, wilcox_p
    )
  )

# --- 4) Build mapping for discrete scales ---
lab_map <- setNames(legend_info$label, legend_info$Level)
legend_breaks <- levels(slope_all$Level)

# --- 5) Plot with updated legend labels (save main + legend separately) ---
suppressPackageStartupMessages({
  library(cowplot)
  library(ggplotify)
})

# Base plot (without legend)
p_main <- ggplot(
  slope_all,
  aes(x = Max_minus_Min_Known_Lifespan, y = slope,
      color = Level, shape = Level)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "black") +
  geom_point(size = 3, alpha = 0.9) +
  theme_all_black +
  theme(
    axis.title.x  = element_text(size = 30, color = "black"),
    axis.title.y  = element_text(size = 30, color = "black"),
    axis.text.x   = element_text(size = 26, color = "black"),
    axis.text.y   = element_text(size = 26, color = "black"),
    legend.position = "none"
  ) +
  labs(
    x = "Lifespan range within group (ln[max/min])",
    y = "Slope (Predicted ~ Known lifespan)"
  ) +
  scale_color_discrete(breaks = legend_breaks, labels = lab_map[legend_breaks]) +
  scale_shape_discrete(breaks = legend_breaks, labels = lab_map[legend_breaks])

ggsave("scatter_slope_vs_lifespan_range_main.png",
       p_main, width = 12, height = 8, dpi = 300)

# Legend only (with roomier spacing)
p_with_legend <- p_main +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 22, color = "black"),
    legend.text       = element_text(size = 18, color = "black", lineheight = 1.6),
    legend.key.height = unit(26, "pt"),
    legend.key.width  = unit(20, "pt")
  )

legend_only <- cowplot::get_legend(p_with_legend)

ggsave("scatter_slope_vs_lifespan_range_legend.png",
       legend_only, width = 6.8, height = 8.6, dpi = 300, bg = "white")









# slope as response
m1 <- lm(slope ~ Level + Max_minus_Min_Known_Lifespan, data = slope_all)

library(car)
anova_res <- Anova(m1, type = "II")   # gives sums of squares type II
anova_res

# Save as CSV
write.csv(as.data.frame(anova_res),
          file = "anova_results_for_slopes.csv",
          row.names = TRUE)

