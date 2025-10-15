setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_final_model/01_01_Fig1_lifespan_scatter_and_bar_plots")

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

# Colors for Train/Test
COLORS <- c("Train" = "#00BFC4", "Test" = "#F8766D")

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

# Keep a species/name column if present (optional)
species_col <- intersect(c("Organism.Name","Organism_Name","Species","Final_Tip_Name"), names(df))
species_col <- if (length(species_col) > 0) species_col[1] else NULL

# ===================== Compute errors =====================
# Only Train & Test
df_tt <- df %>%
  filter(Dataset %in% c("Train", "Test")) %>%
  mutate(
    # Back-transform from ln(days) -> years
    Known_days  = exp(Known_Log_Lifespan),
    Pred_days   = exp(Predicted_Log_Lifespan_Model_Averaging),
    Known_years = Known_days / 365,
    Pred_years  = Pred_days  / 365,

    # Convert to ln(years)
    Known_ln_years = Known_Log_Lifespan - log(365),
    Pred_ln_years  = Predicted_Log_Lifespan_Model_Averaging - log(365),

    # Absolute Error in years
    Abs_Error_years = abs(Pred_years - Known_years),

    # Relative Error on ln(years) scale (%)
    Log_Diff = Pred_ln_years - Known_ln_years,
    Rel_Error_pct_log_units = ifelse(
      abs(Known_ln_years) < 1e-12,
      NA_real_,
      abs(Log_Diff) / abs(Known_ln_years) * 100
    ),

    # *** Force plotting order: Train left, Test right
    Dataset = factor(Dataset, levels = c("Train", "Test"))
  )


# Columns to write
out_cols <- unique(c(
  "Dataset", species_col,
  "Known_years","Pred_years",
  "Abs_Error_years",
  "Rel_Error_pct_log_units",
  "Known_Log_Lifespan","Predicted_Log_Lifespan_Model_Averaging"
))

errors_out <- df_tt %>%
  dplyr::select(all_of(out_cols)) %>%
  dplyr::arrange(Dataset, desc(Abs_Error_years))

out_csv <- file.path(output_dir, "errors_train_test_by_species.csv")
write_csv(errors_out, out_csv)
message("Wrote: ", out_csv)

# ===================== Summaries =====================
summary_tbl <- df_tt %>%
  group_by(Dataset) %>%
  summarise(
    n = dplyr::n(),
    Mean_Abs_Error_years   = mean(Abs_Error_years, na.rm = TRUE),
    Median_Abs_Error_years = median(Abs_Error_years, na.rm = TRUE),
    Mean_Rel_Error_pct_log_units   = mean(Rel_Error_pct_log_units, na.rm = TRUE),
    Median_Rel_Error_pct_log_units = median(Rel_Error_pct_log_units, na.rm = TRUE),
    .groups = "drop"
  )

summary_csv <- file.path(output_dir, "errors_train_test_summary.csv")
write_csv(summary_tbl, summary_csv)
message("Wrote: ", summary_csv)

# ===================== Plots =====================

# ===================== Metrics by Dataset =====================
# Calculate MSE and R² on ln(years) scale
df_tt <- df_tt %>%
  mutate(
    Known_ln_years = Known_Log_Lifespan - log(365),
    Pred_ln_years  = Predicted_Log_Lifespan_Model_Averaging - log(365)
  )

metrics_tbl <- df_tt %>%
  group_by(Dataset) %>%
  summarise(
    MSE = mean((Pred_ln_years - Known_ln_years)^2, na.rm = TRUE),
    R2  = 1 - sum((Pred_ln_years - Known_ln_years)^2, na.rm = TRUE) /
              sum((Known_ln_years - mean(Known_ln_years, na.rm = TRUE))^2, na.rm = TRUE),
    .groups = "drop"
  )

print(metrics_tbl)



# 1) Known vs Predicted (ln years)
p_scatter_all_years <- ggplot(df_tt,
         aes(x = Known_Log_Lifespan - log(365),
             y = Predicted_Log_Lifespan_Model_Averaging - log(365),
             color = Dataset)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2, linewidth = 0.5, color = "grey40") +
  geom_point(alpha = 0.8, size = 2, shape = 16) +
  scale_color_manual(values = COLORS) +
  scale_x_continuous(breaks = seq(-4, 6, 1), limits = c(-4, 6)) +
  scale_y_continuous(breaks = seq(-4, 6, 1), limits = c(-4, 6)) +
  coord_equal(xlim = c(-4, 6), ylim = c(-4, 6)) +
  labs(
    x = "ln(Known lifespan years)",
    y = "ln(Predicted lifespan years)"
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")


# Create label text for Train and Test
metrics_labels <- metrics_tbl %>%
  mutate(
    label = paste0(Dataset, ": MSE=", round(MSE, 2), ", R²=", round(R2, 3))
  )

# Add annotation at different positions
p_scatter_all_years <- p_scatter_all_years +
  annotate("text", x = -3.8, y = 5.8, hjust = 0, vjust = 1,
           label = metrics_labels$label[metrics_labels$Dataset=="Train"],
           size = 4, color = COLORS["Train"]) +
  annotate("text", x = -3.8, y = 5.1, hjust = 0, vjust = 1,
           label = metrics_labels$label[metrics_labels$Dataset=="Test"],
           size = 4, color = COLORS["Test"])

print(p_scatter_all_years)

save_plot(p_scatter_all_years, "plot_known_vs_pred_ln_years_train_test.png", w = 6, h = 5)



# ===================== Fixed label placement =====================

# Absolute error stats (means only)
means_abs <- df_tt %>%
  group_by(Dataset) %>%
  summarise(mean_abs = mean(Abs_Error_years, na.rm = TRUE), .groups = "drop")

# Relative error stats (means only)
means_rel <- df_tt %>%
  group_by(Dataset) %>%
  summarise(mean_rel = mean(Rel_Error_pct_log_units, na.rm = TRUE), .groups = "drop")

# 2) Absolute Error plot (violin version)
p_violin_abs_log <- ggplot(df_tt, aes(x = Dataset, y = pmax(Abs_Error_years, 1e-6), fill = Dataset)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  geom_text(
    data = means_abs,
    aes(x = Dataset,
        y = 120,  # fixed position at y = 120
        label = paste0("Mean = ", sprintf("%.2f", mean_abs), " y"),
        color = Dataset),
    inherit.aes = FALSE, size = 5, fontface = "bold"
  ) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c("0.01", "0.1", "1", "10", "100"),
                limits = c(0.01, 120)) +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  labs(x = "", y = "Absolute error (years)") +
  theme_minimal(base_size = 19) +
  theme(legend.position = "none")

save_plot(p_violin_abs_log, "violin_abs_error_years_log10_train_test.png", w = 5, h = 5)



# 3) Relative Error plot (violin version)
p_violin_rel_log_units <- ggplot(df_tt, aes(x = Dataset, y = Rel_Error_pct_log_units, fill = Dataset)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8, na.rm = TRUE) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1, na.rm = TRUE) +
  geom_text(
    data = means_rel,
    aes(x = Dataset,
        y = 1200,
        label = paste0("Mean = ", sprintf("%.1f", mean_rel), "%"),
        color = Dataset),
    inherit.aes = FALSE, size = 5, fontface = "bold"
  ) +
  scale_y_continuous(breaks = seq(0, 1200, 200), limits = c(0, 1200)) +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  labs(x = "", y = "Relative error (%)") +
  theme_minimal(base_size = 19) +
  theme(legend.position = "none")

save_plot(p_violin_rel_log_units, "violin_rel_error_pct_ln_years_train_test.png", w = 5, h = 5)



# ===================== Console outputs =====================
print(summary_tbl)



