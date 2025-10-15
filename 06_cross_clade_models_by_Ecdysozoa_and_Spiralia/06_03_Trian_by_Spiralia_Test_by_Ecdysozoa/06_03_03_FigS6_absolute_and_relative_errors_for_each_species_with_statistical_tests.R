setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/06_cross_clade_models_by_Ecdysozoa_and_Spiralia/06_03_Trian_by_Spiralia_Test_by_Ecdysozoa")

# ===================== Setup =====================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(ggpubr) 
})

input_path <- "log_lifespan_predictions.csv"
output_dir <- dirname(input_path)

# Colors for groups
COLORS <- c("Spiralia (Train)" = "#6a3d9a",
            "Ecdysozoa (Test)" = "#33a02c")

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
df_tt <- df %>%
  filter(Dataset %in% c("Train", "Test")) %>%
  mutate(
    # Rename datasets
    Dataset = dplyr::recode(Dataset,
                     "Train" = "Spiralia (Train)",
                     "Test"  = "Ecdysozoa (Test)"),

    # Back-transform from ln(days) -> years
    Known_days  = exp(Known_Log_Lifespan),
    Pred_days   = exp(Predicted_Log_Lifespan_Model_Averaging),
    Known_years = Known_days / 365,
    Pred_years  = Pred_days  / 365,

    # Absolute Error in years
    Abs_Error_years = abs(Pred_years - Known_years),

    # Relative Error on ln-days scale (%)
    Log_Diff = Predicted_Log_Lifespan_Model_Averaging - Known_Log_Lifespan,
    Rel_Error_pct_log_units = ifelse(
      abs(Known_Log_Lifespan) < 1e-12,
      NA_real_,
      abs(Log_Diff) / abs(Known_Log_Lifespan) * 100
    ),

    # Force plotting order
    Dataset = factor(Dataset, levels = c("Spiralia (Train)", "Ecdysozoa (Test)"))
  )

# ===================== Write errors =====================
out_cols <- unique(c(
  "Dataset", species_col,
  "Known_years","Pred_years",
  "Abs_Error_years",
  "Rel_Error_pct_log_units",
  "Known_Log_Lifespan","Predicted_Log_Lifespan_Model_Averaging"
))

errors_out <- df_tt %>%
  select(all_of(out_cols)) %>%
  arrange(Dataset, desc(Abs_Error_years))

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

# ===================== Metrics by Dataset =====================
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

# ===================== Scatter plot =====================
p_scatter_all_years <- ggplot(df_tt,
         aes(x = Known_ln_years,
             y = Pred_ln_years,
             color = Dataset)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2, linewidth = 0.5, color = "grey40") +
  geom_point(alpha = 0.8, size = 2, shape = 16) +
  scale_color_manual(values = COLORS) +
  scale_x_continuous(breaks = seq(-4, 6, 1), limits = c(-4, 6)) +
  scale_y_continuous(breaks = seq(-4, 6, 1), limits = c(-4, 6)) +
  coord_equal(xlim = c(-4, 6), ylim = c(-4, 6)) +
  labs(x = "ln(Known lifespan years)",
       y = "ln(Predicted lifespan years)") +
  theme_minimal(base_size = 17) +
  theme(legend.position = "none")

metrics_labels <- metrics_tbl %>%
  mutate(label = paste0(Dataset, ": MSE=", round(MSE, 2), ", R²=", round(R2, 3)))

p_scatter_all_years <- p_scatter_all_years +
  annotate("text", x = -3.8, y = 5.8, hjust = 0, vjust = 1,
           label = metrics_labels$label[metrics_labels$Dataset=="Spiralia (Train)"],
           size = 4, color = COLORS["Spiralia (Train)"]) +
  annotate("text", x = -3.8, y = 5.1, hjust = 0, vjust = 1,
           label = metrics_labels$label[metrics_labels$Dataset=="Ecdysozoa (Test)"],
           size = 4, color = COLORS["Ecdysozoa (Test)"])

save_plot(p_scatter_all_years, "plot_known_vs_pred_ln_years_train_test.png", w = 6, h = 5)
print(p_scatter_all_years)


# ===================== Correlations (Pearson's r + global FDR) =====================
library(stats)

# Per-dataset Pearson r and p (remove NA pairs first)
cor_tbl <- df_tt %>%
  group_by(Dataset) %>%
  summarise(
    r = cor(Known_ln_years, Pred_ln_years, use = "pairwise.complete.obs", method = "pearson"),
    p = {
      tmp <- na.omit(data.frame(x = Known_ln_years, y = Pred_ln_years))
      if (nrow(tmp) >= 3) cor.test(tmp$x, tmp$y, method = "pearson")$p.value else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(FDR_p = p.adjust(p, method = "fdr"),
         label = paste0("r = ", round(r, 3),
                        ", FDR p = ", format(FDR_p, scientific = TRUE, digits = 3)))

print(cor_tbl)

# ===================== Add labels to the scatter plot =====================
p_scatter_all_years <- p_scatter_all_years +
  annotate("text", x = -3.8, y = 4.4, hjust = 0, vjust = 1,
           label = cor_tbl$label[cor_tbl$Dataset == "Spiralia (Train)"],
           size = 4, color = COLORS["Spiralia (Train)"]) +
  annotate("text", x = -3.8, y = 3.7, hjust = 0, vjust = 1,
           label = cor_tbl$label[cor_tbl$Dataset == "Ecdysozoa (Test)"],
           size = 4, color = COLORS["Ecdysozoa (Test)"])

save_plot(p_scatter_all_years, "plot_known_vs_pred_ln_years_with_r_FDR.png", w = 6, h = 5)
print(p_scatter_all_years)


# ===================== Wilcoxon Tests =====================
p_known <- wilcox.test(Known_ln_years ~ Dataset, data = df_tt)$p.value
p_pred  <- wilcox.test(Pred_ln_years  ~ Dataset, data = df_tt)$p.value
p_rel   <- wilcox.test(Rel_Error_pct_log_units ~ Dataset,
                       data = df_tt %>% filter(!is.na(Rel_Error_pct_log_units)))$p.value

message("Wilcoxon p (Known ln-years): ", format(p_known, scientific = TRUE, digits = 3))
message("Wilcoxon p (Pred  ln-years): ", format(p_pred,  scientific = TRUE, digits = 3))
message("Wilcoxon p (Relative error %): ", format(p_rel, scientific = TRUE, digits = 3))

p_known_label <- paste0("p = ", format(p_known, scientific = TRUE, digits = 3))
p_pred_label  <- paste0("p = ", format(p_pred,  scientific = TRUE, digits = 3))
p_rel_label   <- paste0("p = ", format(p_rel,   scientific = TRUE, digits = 3))

# ===================== Box plots =====================
# Predicted lifespan
p_box_pred_ln_years <- ggplot(df_tt, aes(x = Dataset, y = Pred_ln_years, fill = Dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  scale_fill_manual(values = COLORS) +
  scale_x_discrete(labels = c("Spiralia (Train)" = "Spiralia\n(Train)",
                              "Ecdysozoa (Test)" = "Ecdysozoa\n(Test)")) +
  scale_y_continuous(breaks = seq(-4, 6, 1), limits = c(-4, 6.5)) +
  labs(x = "", y = "ln(Predicted lifespan years)") +
  theme_minimal(base_size = 17) +
  theme(legend.position = "none") +
  annotate("segment", x = 1, xend = 2, y = 5.9, yend = 5.9, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = 5.9, yend = 5.7, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = 5.9, yend = 5.7, linewidth = 0.4) +
  annotate("text", x = 1.5, y = 6.4, label = p_pred_label, size = 5)

save_plot(p_box_pred_ln_years, "box_predicted_lifespan_ln_years.png", w = 3, h = 5)
print(p_box_pred_ln_years)


# Known lifespan
p_box_known_ln_years <- ggplot(df_tt, aes(x = Dataset, y = Known_ln_years)) +
  geom_boxplot(fill = "grey80", color = "black", outlier.shape = NA, width = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  scale_fill_manual(values = COLORS) +
  scale_x_discrete(labels = c("Spiralia (Train)" = "Spiralia",
                              "Ecdysozoa (Test)" = "Ecdysozoa")) +
  scale_y_continuous(breaks = seq(-4, 6, 1), limits = c(-4, 6.5)) +
  labs(x = "", y = "ln(Known lifespan years)") +
  theme_minimal(base_size = 17) +
  theme(legend.position = "none") +
  annotate("segment", x = 1, xend = 2, y = 5.9, yend = 5.9, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = 5.9, yend = 5.7, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = 5.9, yend = 5.7, linewidth = 0.4) +
  annotate("text", x = 1.5, y = 6.4, label = p_known_label, size = 5)

save_plot(p_box_known_ln_years, "box_known_lifespan_ln_years_greyscale.png", w = 3, h = 5)
print(p_box_known_ln_years)


# Relative error
p_box_rel_log_units <- ggplot(df_tt, aes(x = Dataset, y = Rel_Error_pct_log_units, fill = Dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, na.rm = TRUE) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1, na.rm = TRUE) +
  scale_fill_manual(values = COLORS) +
  scale_x_discrete(labels = c("Spiralia (Train)" = "Spiralia\n(Train)",
                              "Ecdysozoa (Test)" = "Ecdysozoa\n(Test)")) +
  scale_y_continuous(breaks = seq(0, 120, 20), limits = c(0, 120)) +
  labs(x = "", y = "Relative error (%)") +
  theme_minimal(base_size = 17) +
  theme(legend.position = "none") +
  annotate("segment", x = 1, xend = 2, y = 112, yend = 112, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = 112, yend = 110, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = 112, yend = 110, linewidth = 0.4) +
  annotate("text", x = 1.5, y = 119, label = p_rel_label, size = 5)

save_plot(p_box_rel_log_units, "box_rel_error_pct_ln_days.png", w = 3, h = 5)
print(p_box_rel_log_units)


# ===================== Console outputs =====================
print(summary_tbl)
