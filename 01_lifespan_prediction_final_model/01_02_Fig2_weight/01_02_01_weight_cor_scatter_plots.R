# ======================= Setup =======================
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_and_error/01_02_Fig2_weight")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
})

# ----------------------- Config ----------------------
LEGEND_FILE   <- "legend_region_caption_vertical.png"
PLOT1_FILE    <- "scatter_Correlation_vs_AverageWeight.png"
PLOT2_FILE    <- "scatter_AvgWeightPerFold_vs_SelectionFrequency.png"
PLOT3_FILE    <- "scatter_Correlation_vs_SelectionFrequency.png"

FIG_W <- 8
FIG_H <- 6
DPI   <- 300

LABEL_FREQ_THRESH <- 0.5   # label features with Selection_Frequency >= this

REGION_LEVELS <- c("Genome", "Upstream", "Exon", "Intron", "Downstream")

region_colors <- c(
  "Genome"     = "#1b9e77",
  "Upstream"   = "#d95f02",
  "Exon"       = "#7570b3",
  "Intron"     = "#e7298a",
  "Downstream" = "#66a61e"
)

# ------------------ Helper functions -----------------
parse_region <- function(feature_name) {
  raw <- str_extract(feature_name, "(?<=_Rate_).*")
  raw_lower <- tolower(raw)
  out <- dplyr::case_when(
    str_detect(raw_lower, "upstream")   ~ "Upstream",
    str_detect(raw_lower, "downstream") ~ "Downstream",
    str_detect(raw_lower, "intron")     ~ "Intron",
    str_detect(raw_lower, "exon")       ~ "Exon",
    str_detect(raw_lower, "genome")     ~ "Genome",
    TRUE ~ NA_character_
  )
  factor(out, levels = REGION_LEVELS)
}

make_label_tetranuc <- function(feature_name) {
  str_extract(feature_name, "^[ACGT]+")
}

compute_avg_weight_per_selected_fold <- function(df, fold_cols) {
  df %>%
    rowwise() %>%
    mutate(
      Sum_Weight    = sum(c_across(all_of(fold_cols)), na.rm = TRUE),
      Nonzero_Folds = sum(c_across(all_of(fold_cols)) != 0, na.rm = TRUE),
      Avg_Weight_Per_Fold = ifelse(Nonzero_Folds > 0, Sum_Weight / Nonzero_Folds, 0)
    ) %>%
    ungroup()
}

# ======================= Load ========================
df_raw <- read_csv("combined_model_weights.csv", show_col_types = FALSE)

# Ensure FDR column exists
stopifnot("FDR_adjusted_p" %in% names(df_raw))

# Common enrichment (used by all plots)
df_enriched <- df_raw %>%
  filter(Feature_Name != "(Intercept)") %>%
  mutate(
    Region_short = parse_region(Feature_Name),
    Tetranuc     = make_label_tetranuc(Feature_Name),
    Label        = Tetranuc,
    Significant  = !is.na(FDR_adjusted_p) & FDR_adjusted_p <= 0.05
  ) %>%
  filter(!is.na(Region_short))  # ensure only the 5 regions remain

# Fold columns (if present)
fold_cols <- grep("^Fold_", names(df_enriched), value = TRUE)

# ======================= Plot 1 =======================
# Correlation (x) vs Average_Weight (y), colored by region, alpha by significance
p1 <- ggplot(df_enriched,
             aes(x = Lifespan_Feature_Correlation_All_Species,
                 y = Average_Weight,
                 color = Region_short,
                 alpha = Significant)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = region_colors, drop = TRUE) +
  scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.15), guide = "none") +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  labs(x = "Pearson's correlation (r)", y = "Average weight")

ggsave(PLOT1_FILE, p1, width = FIG_W, height = FIG_H, dpi = DPI)

# ======================= Plot 2 =======================
# Avg_Weight_Per_Fold (x) vs Selection_Frequency (y) with labels >= threshold
if (length(fold_cols) > 0) {
  df_avgfold <- compute_avg_weight_per_selected_fold(df_enriched, fold_cols) %>%
    mutate(Region_short = factor(Region_short, levels = REGION_LEVELS))

  p2 <- ggplot(df_avgfold,
               aes(x = Avg_Weight_Per_Fold,
                   y = Selection_Frequency,
                   color = Region_short)) +
    geom_point(alpha = 0.7, size = 2.3) +
    geom_text_repel(
      data = subset(df_avgfold, Selection_Frequency >= LABEL_FREQ_THRESH),
      aes(label = Label, color = Region_short),
      size = 5,
      max.overlaps = Inf,
      box.padding = 0.6,
      point.padding = 0.4,
      force = 2,
      min.segment.length = 0,
      show.legend = FALSE
    ) +
    scale_color_manual(values = region_colors, drop = TRUE) +
    theme_minimal(base_size = 22) +
    theme(legend.position = "none") +
    labs(x = "Average weight per selected fold", y = "Selection frequency")

  ggsave(PLOT2_FILE, p2, width = FIG_W, height = FIG_H, dpi = DPI)

  # --------- Save vertical shared legend (ordered) ---------
  legend_grob <- get_legend(
    ggplot(df_avgfold, aes(x = Avg_Weight_Per_Fold,
                           y = Selection_Frequency,
                           color = Region_short)) +
      geom_point(size = 4) +
      scale_color_manual(values = region_colors, drop = TRUE) +
      theme_minimal(base_size = 28) +
      theme(
        legend.position = "right",        # vertical
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 24)
      ) +
      labs(color = "Region")
  )
  ggsave(LEGEND_FILE, ggdraw(legend_grob),
         width = 4, height = 8, dpi = DPI)
} else {
  message("No fold columns detected (starting with 'Fold_'); skipping Plot 2 and legend.")
}

# ======================= Plot 3 =======================
# Correlation (x) vs Selection_Frequency (y), colored by region, alpha by significance
p3 <- ggplot(df_enriched,
             aes(x = Lifespan_Feature_Correlation_All_Species,
                 y = Selection_Frequency,
                 color = Region_short,
                 alpha = Significant)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = region_colors, drop = TRUE) +
  scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.15), guide = "none") +
  theme_minimal(base_size = 22) +
  theme(legend.position = "none") +
  labs(x = "Pearson's correlation (r)", y = "Selection frequency")

ggsave(PLOT3_FILE, p3, width = FIG_W, height = FIG_H, dpi = DPI)

# ======================= Done ========================
message("Wrote: ", PLOT1_FILE)
if (length(fold_cols) > 0) {
  message("Wrote: ", PLOT2_FILE)
  message("Wrote: ", LEGEND_FILE)
}
message("Wrote: ", PLOT3_FILE)
