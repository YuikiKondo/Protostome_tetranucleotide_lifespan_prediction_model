setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_and_error/01_05_FigS1_absolute_and_relative_errors_by_knwon_lifespan_group")

# ===================== Setup =====================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
})

input_path <- "log_lifespan_predictions.csv"
output_dir <- dirname(input_path)

# Colors (Train = orange, Test = blue-green)
COLORS <- c("Train" = "#F8766D", "Test" = "#00BFC4")

# Fixed y-positions for mean labels
y_signed_label_close <- 15
y_signed_label_broad <- 15
y_abs_label          <- 120
y_rel_label          <- 120

# Helper: safe saver
save_plot <- function(p, fname, w = 7, h = 5, dpi = 300) {
  out <- file.path(output_dir, fname)
  ggsave(filename = out, plot = p, width = w, height = h, dpi = dpi)
  message("Wrote: ", out)
}

# ===================== Load =====================
df <- read_csv(input_path, show_col_types = FALSE)

req_cols <- c("Dataset", "Known_Log_Lifespan", "Predicted_Log_Lifespan_Model_Averaging")
missing <- setdiff(req_cols, names(df))
if (length(missing) > 0) stop("Missing required column(s): ", paste(missing, collapse = ", "))

# ===================== Compute errors =====================
df_tt <- df %>%
  filter(Dataset %in% c("Train", "Test")) %>%
  mutate(
    # back-transform
    Known_days  = exp(Known_Log_Lifespan),
    Pred_days   = exp(Predicted_Log_Lifespan_Model_Averaging),
    Known_years = Known_days / 365,
    Pred_years  = Pred_days  / 365,

    # errors
    Error_years      = Pred_years - Known_years,
    Abs_Error_years  = abs(Error_years),
    Log_Diff         = Predicted_Log_Lifespan_Model_Averaging - Known_Log_Lifespan,
    Rel_Error_pct_log_units = abs(Log_Diff) / abs(Known_Log_Lifespan) * 100,

    Dataset = factor(Dataset, levels = c("Train", "Test")),

    # ===================== Lifespan bins by YEARS =====================
    Lifespan_years_bin = cut(
      Known_years,
      breaks = c(-Inf, 0.5, 2, 5, 20, 40, Inf),
      labels = c("<0.5", "0.5–<2", "2–<5", "5–<20", "20–<40", "≥40"),
      right = FALSE
    )
  )

# ===================== Group means for labels =====================
means_df <- df_tt %>%
  filter(!is.na(Lifespan_years_bin)) %>%
  group_by(Lifespan_years_bin, Dataset) %>%
  summarise(
    mean_signed = mean(Error_years, na.rm = TRUE),
    mean_abs    = mean(Abs_Error_years, na.rm = TRUE),
    mean_rel    = mean(Rel_Error_pct_log_units, na.rm = TRUE),
    .groups = "drop"
  )

# A little headroom so top labels don't clip
top_margin_theme <- theme(
  plot.margin = margin(t = 16, r = 10, b = 10, l = 10),
  legend.position = "right"
)

# ===================== Signed error (close-up) =====================
p_error_bin_box_close <- ggplot(df_tt %>% filter(!is.na(Lifespan_years_bin)),
  aes(x = Lifespan_years_bin, y = Error_years, fill = Dataset)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_boxplot(outlier.shape = NA, width = 0.7,
               position = position_dodge2(width = 0.9)) +
  geom_jitter(aes(group = Dataset),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.35, size = 1) +
  geom_text(
    data = means_df,
    aes(x = Lifespan_years_bin, y = y_signed_label_close,
        label = sprintf("%.2f", mean_signed), color = Dataset),
    position = position_dodge2(width = 1), vjust = -0.25, size = 4
  ) +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  scale_y_continuous(limits = c(-10, 15), expand = expansion(mult = c(0, 0.10))) +
  labs(x = "Known lifespan (years)", y = "Error (years)") +
  theme_minimal(base_size = 16) + top_margin_theme

save_plot(p_error_bin_box_close, "box_error_years_by_year_bins_close.png", w = 8, h = 5.2)

# ===================== Signed error (broad) =====================
p_error_bin_box_broad <- ggplot(df_tt %>% filter(!is.na(Lifespan_years_bin)),
  aes(x = Lifespan_years_bin, y = Error_years, fill = Dataset)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_boxplot(outlier.shape = NA, width = 0.7,
               position = position_dodge2(width = 0.9)) +
  geom_jitter(aes(group = Dataset),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.35, size = 1) +
  geom_text(
    data = means_df,
    aes(x = Lifespan_years_bin, y = y_signed_label_broad,
        label = sprintf("%.2f", mean_signed), color = Dataset),
    position = position_dodge2(width = 1), vjust = -0.25, size = 4
  ) +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  scale_y_continuous(limits = c(-250, 15), expand = expansion(mult = c(0, 0.10))) +
  labs(x = "Known lifespan (years)", y = "Error (years)") +
  theme_minimal(base_size = 16) + top_margin_theme

save_plot(p_error_bin_box_broad, "box_error_years_by_year_bins_broad.png", w = 8, h = 5.2)

# ===================== Absolute error (years, log10) =====================
p_abs_bin_box <- ggplot(df_tt %>% filter(!is.na(Lifespan_years_bin)),
  aes(x = Lifespan_years_bin, y = pmax(Abs_Error_years, 1e-6), fill = Dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.7,
               position = position_dodge2(width = 0.9)) +
  geom_jitter(aes(group = Dataset),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.35, size = 1) +
  geom_text(
    data = means_df,
    aes(x = Lifespan_years_bin, y = y_abs_label,
        label = sprintf("%.2f", mean_abs), color = Dataset),
    position = position_dodge2(width = 1), vjust = -0.25, size = 4
  ) +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01","0.1","1","10","100"),
    expand = expansion(mult = c(0.02, 0.25))
  ) +
  labs(x = "Known lifespan (years)", y = "Absolute error (years)") +
  theme_minimal(base_size = 16) + top_margin_theme

save_plot(p_abs_bin_box, "box_abs_error_years_by_year_bins_log10.png", w = 8, h = 5.2)

# ===================== Relative error (%) on ln-days, log10 =====================
pct_breaks <- c(0.01, 0.1, 1, 10, 100)
pct_labels <- paste0(pct_breaks, "%")

p_rel_bin_box <- ggplot(df_tt %>% filter(!is.na(Lifespan_years_bin)),
  aes(x = Lifespan_years_bin, y = pmax(Rel_Error_pct_log_units, 1e-4), fill = Dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.7,
               position = position_dodge2(width = 0.9), na.rm = TRUE) +
  geom_jitter(aes(group = Dataset),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.35, size = 1, na.rm = TRUE) +
  geom_text(
    data = means_df,
    aes(x = Lifespan_years_bin, y = y_rel_label,
        label = sprintf("%.2f", mean_rel), color = Dataset),
    position = position_dodge2(width = 1), vjust = -0.25, size = 4
  ) +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  scale_y_log10(
    breaks = pct_breaks, labels = pct_labels,
    expand = expansion(mult = c(0.02, 0.25))
  ) +
  labs(x = "Known lifespan (years)", y = "Relative error (%)") +
  theme_minimal(base_size = 16) + top_margin_theme

save_plot(p_rel_bin_box, "box_rel_error_pct_ln_days_by_year_bins_log10.png", w = 8, h = 5.2)
