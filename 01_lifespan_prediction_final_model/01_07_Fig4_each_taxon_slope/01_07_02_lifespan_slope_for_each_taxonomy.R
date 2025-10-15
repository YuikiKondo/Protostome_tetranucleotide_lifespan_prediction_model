# ====== Set working directory ======
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_final_model/01_07_Fig4_each_taxon_slope")

# ====== Load required libraries ======
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(broom)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# ---------- Reusable all-black (not bold) theme ----------
theme_all_black <- theme_classic(base_size = 18) +
  theme(
    plot.title        = element_text(color = "black"),
    axis.title        = element_text(color = "black"),
    axis.text.x       = element_text(color = "black", size = 18),
    axis.text.y       = element_text(color = "black", size = 18),
    legend.title      = element_text(color = "black"),
    legend.text       = element_text(color = "black"),
    strip.text        = element_text(color = "black"),
    axis.line         = element_line(color = "black", linewidth = 0.6),
    axis.ticks        = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(4, "pt"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )

# ====== Step 1: Read data ======
predictions <- read.csv("log_lifespan_predictions_with_taxonomy.csv")

# ====== Step 1.5: Merge Arthropoda groups ======
predictions <- predictions %>%
  mutate(
    Phylum = case_when(
      Phylum %in% c("Arthropoda (Chelicerata)", "Arthropoda (Mandibulata)") ~ "Arthropoda",
      TRUE ~ Phylum
    )
  )

# ====== NEW: switch to ln(years) from ln(days) ======
# ln(years) = ln(days) - ln(365)
L365 <- log(365)
predictions <- predictions %>%
  mutate(
    Known_Log_Lifespan_years     = Known_Log_Lifespan - L365,
    Predicted_Log_Lifespan_years = Predicted_Log_Lifespan_Model_Averaging - L365
  )

# ====== Step 2: Filter for Train and Test datasets ======
filtered_data <- predictions %>%
  filter(Dataset %in% c("Train", "Test"))

# ====== Step 3: Define taxonomic levels ======
taxonomic_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

# ====== Step 4: Define function to compute slope per group (using ln-years) ======
compute_slope_by_level <- function(level) {
  higher_levels <- taxonomic_levels[1:which(taxonomic_levels == level) - 1]

  # 1) Canonical (most frequent) lineage per focal group
  if (length(higher_levels) > 0) {
    lineage_table <- filtered_data %>%
      filter(!is.na(.data[[level]])) %>%
      count(across(all_of(c(level, higher_levels))), name = "n") %>%
      arrange(desc(n)) %>%
      group_by(across(all_of(level))) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(-n)
  } else {
    lineage_table <- filtered_data %>%
      filter(!is.na(.data[[level]])) %>%
      distinct(across(all_of(level)))
  }

  # 2) Slope per focal group (ln-years)
  result <- filtered_data %>%
    filter(!is.na(.data[[level]])) %>%
    group_by(across(all_of(level))) %>%
    mutate(n_species = n()) %>%
    ungroup() %>%
    filter(n_species > 2) %>%
    group_by(across(all_of(level))) %>%
    nest() %>%
    mutate(
      n_species = map_int(data, nrow),
      model = map(data, ~ lm(Predicted_Log_Lifespan_years ~ Known_Log_Lifespan_years, data = .x)),
      slope = map_dbl(model, ~ coef(.x)[["Known_Log_Lifespan_years"]]),
      mse = map_dbl(data, ~ mean((.x$Predicted_Log_Lifespan_years - .x$Known_Log_Lifespan_years)^2)),
      mae = map_dbl(data, ~ mean(abs(.x$Predicted_Log_Lifespan_years - .x$Known_Log_Lifespan_years))),
      mean_error = map_dbl(data, ~ mean(.x$Predicted_Log_Lifespan_years - .x$Known_Log_Lifespan_years)),
      Max_Known_Lifespan = map_dbl(data, ~ max(.x$Known_Log_Lifespan_years, na.rm = TRUE)),
      Min_Known_Lifespan = map_dbl(data, ~ min(.x$Known_Log_Lifespan_years, na.rm = TRUE)),
      Max_minus_Min_Known_Lifespan = Max_Known_Lifespan - Min_Known_Lifespan
    ) %>%
    select(all_of(level), n_species, slope, mse, mae, mean_error,
           Max_Known_Lifespan, Min_Known_Lifespan, Max_minus_Min_Known_Lifespan) %>%
    ungroup()

  # 3) Join single canonical lineage
  result <- result %>%
    left_join(lineage_table, by = level) %>%
    mutate(Level = level) %>%
    relocate(all_of(c(level, higher_levels)), .before = n_species)

  result %>% distinct(across(all_of(level)), .keep_all = TRUE)
}

# ====== Step 4.5: Compute slope for ALL species ("Protostomia") in ln-years ======
compute_slope_all <- function() {
  dat <- filtered_data %>%
    select(Known_Log_Lifespan_years, Predicted_Log_Lifespan_years)

  n_species <- nrow(dat)
  mdl <- lm(Predicted_Log_Lifespan_years ~ Known_Log_Lifespan_years, data = dat)

  tibble::tibble(
    Level = "Protostomia",
    n_species = n_species,
    slope = coef(mdl)[["Known_Log_Lifespan_years"]],
    mse = mean((dat$Predicted_Log_Lifespan_years - dat$Known_Log_Lifespan_years)^2),
    mae = mean(abs(dat$Predicted_Log_Lifespan_years - dat$Known_Log_Lifespan_years)),
    mean_error = mean(dat$Predicted_Log_Lifespan_years - dat$Known_Log_Lifespan_years),
    Max_Known_Lifespan = max(dat$Known_Log_Lifespan_years, na.rm = TRUE),
    Min_Known_Lifespan = min(dat$Known_Log_Lifespan_years, na.rm = TRUE),
    Max_minus_Min_Known_Lifespan = Max_Known_Lifespan - Min_Known_Lifespan
  )
}

# ====== Step 5: Compute slopes and save CSVs ======
slope_summary_all <- list()

for (level in taxonomic_levels) {
  result <- compute_slope_by_level(level)
  result$Level <- level
  filename <- paste0("slope_by_", tolower(level), "_ln_years.csv")  # renamed to reflect units
  write.csv(dplyr::select(result, -any_of("Level")), filename, row.names = FALSE)
  cat("Saved:", filename, "\n")
  slope_summary_all[[level]] <- result
}

all_res <- compute_slope_all()
write.csv(dplyr::select(all_res, -any_of("Level")),
          "slope_all_protostomia_ln_years.csv", row.names = FALSE)
cat("Saved: slope_all_protostomia_ln_years.csv\n")
slope_summary_all[["Protostomia"]] <- all_res

# ====== Step 6: Combine all into a single dataframe ======
slope_summary_all_df <- bind_rows(slope_summary_all)

# ====== Step 7: Add slope direction ======
slope_summary_all_df <- slope_summary_all_df %>%
  mutate(SlopeDirection = ifelse(is.na(slope), NA, ifelse(slope >= 0, "Positive", "Negative")))

# ====== Step 8: Calculate Positive Proportions (include Protostomia) ======
level_order <- c("Protostomia", taxonomic_levels)

positive_proportion <- slope_summary_all_df %>%
  group_by(Level) %>%
  summarise(
    total_groups = n(),
    positive_groups = sum(SlopeDirection == "Positive", na.rm = TRUE),
    Proportion_Positive = positive_groups / total_groups,
    .groups = "drop"
  ) %>%
  mutate(Level = factor(Level, levels = level_order)) %>%
  arrange(Level)

# ====== Step 9: Save Positive Proportions ======
cat("\nProportion of Positive Slopes by Level:\n")
print(positive_proportion)
write.csv(positive_proportion, "positive_slope_proportions_by_level_ln_years.csv", row.names = FALSE)
cat("Saved: positive_slope_proportions_by_level_ln_years.csv\n")

# ====== Step 10: Wilcoxon tests only (one-sided: greater than zero) ======
suppressWarnings({
  slope_wilcoxon_tests <- slope_summary_all_df %>%
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
          median_slope = stats::median(x, na.rm = TRUE),
          mean_slope = mean(x, na.rm = TRUE),
          mean_mse = mean(.x$mse, na.rm = TRUE),
          wilcox_W = ifelse(is.null(wt), NA_real_, unname(wt$statistic)),
          wilcox_p_value = ifelse(is.null(wt), NA_real_, wt$p.value),
          test_description = "Wilcoxon signed-rank test (median slope > 0, one-sided)"
        )
      } else {
        tibble::tibble(
          n = n,
          median_slope = stats::median(x, na.rm = TRUE),
          mean_slope = mean(x, na.rm = TRUE),
          mean_mse = mean(.x$mse, na.rm = TRUE),
          wilcox_W = NA_real_,
          wilcox_p_value = NA_real_,
          test_description = "Wilcoxon signed-rank test (median slope > 0, one-sided)"
        )
      }
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Level = factor(Level, levels = c("Protostomia", taxonomic_levels))) %>%
    dplyr::arrange(Level)
})

print(slope_wilcoxon_tests)
write.csv(slope_wilcoxon_tests, "slope_wilcoxon_tests_by_level_ln_years.csv", row.names = FALSE)
cat("Saved: slope_wilcoxon_tests_by_level_ln_years.csv\n")

# ====== Step 11: Prepare slope count data for plotting ======
slope_count_plot_data <- slope_summary_all_df %>%
  filter(SlopeDirection %in% c("Positive", "Negative")) %>%
  group_by(Level, SlopeDirection) %>%
  summarise(n = n(), .groups = "drop")

# ====== Step 12: Plot Positive vs Negative Slope Counts (ln-years label) ======
ggplot(slope_count_plot_data, aes(x = Level, y = n, fill = SlopeDirection)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.9),
           width = 0.8,
           color = "black") +
  scale_fill_manual(values = c("Positive" = "#D55E00", "Negative" = "#009E73")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Taxonomic Level", y = "Number of Groups", fill = "Slope Direction") +
  theme_all_black +
  theme(axis.line.x = element_blank())
ggsave("slope_direction_by_taxonomic_level_ln_years.png", width = 10, height = 6)

# ====== Step 13: Plot Proportion of Positive Slopes ======
ggplot(positive_proportion, aes(x = Level, y = Proportion_Positive)) +
  geom_bar(stat = "identity", fill = "grey60", color = "black", width = 0.7) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Taxonomic level", y = "Proportion of positive slope groups") +
  theme_all_black +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text  = element_text(size = 18)
  )


ggsave("positive_slope_proportion_for_each_taxonomic_level_ln_years.png", width = 10, height = 6)



# ===================== Scatter plots (ln-years) =====================

suppressPackageStartupMessages({
  library(cowplot)
  library(ggplotify)
})

extract_legend <- function(p, ncol = 2, inner_left_pt = 28, outer_left_pt = 28) {
  lg <- cowplot::get_legend(
    p +
      guides(color = guide_legend(ncol = ncol, byrow = TRUE)) +
      theme(
        legend.position   = "right",
        legend.box.margin = margin(6, 12, 6, inner_left_pt),
        legend.margin     = margin(0, 0, 0, 0),
        legend.key.width  = unit(14, "pt"),
        legend.key.height = unit(10, "pt"),
        legend.text       = element_text(size = 14, color = "black"),
        legend.title      = element_text(size = 16, color = "black")
      )
  )
  ggplotify::as.ggplot(lg) + theme(plot.margin = margin(4, 12, 4, outer_left_pt))
}

# Helper: nice axis limits based on ln-years data (keeps square plot)
xy_limits_from <- function(df, xcol, ycol, pad = 0.15) {
  x_rng <- range(df[[xcol]], na.rm = TRUE)
  y_rng <- range(df[[ycol]], na.rm = TRUE)
  lo <- floor(min(x_rng[1], y_rng[1]))
  hi <- ceiling(max(x_rng[2], y_rng[2]))
  span <- hi - lo
  lo <- lo - pad * span
  hi <- hi + pad * span
  list(xlim = c(lo, hi), ylim = c(lo, hi))
}

# ====== Axis limits (ln-years) ======
X_MIN <- -4; X_MAX <- 6
Y_MIN <- -4; Y_MAX <- 6
BRKS  <- seq(-4, 6, by = 1)

# ====== Scatter with separate legend (ln-years) ======
plot_scatter_colored_by_group <- function(level,
                                          min_n = 3,
                                          max_groups = Inf,
                                          out_width = 14,
                                          out_height = 10,
                                          dpi = 300) {
  stopifnot(level %in% names(filtered_data))

  df0 <- filtered_data %>%
    filter(!is.na(.data[[level]]))

  group_sizes <- df0 %>%
    count(across(all_of(level)), name = "n_group") %>%
    arrange(desc(n_group)) %>%
    filter(n_group >= min_n)

  if (nrow(group_sizes) == 0) {
    message("No groups with enough points at level: ", level)
    return(invisible(NULL))
  }

  if (is.finite(max_groups) && nrow(group_sizes) > max_groups) {
    group_sizes <- group_sizes %>% slice_head(n = max_groups)
  }
  keep_groups <- group_sizes %>% pull(!!sym(level))

  df_level <- df0 %>%
    filter(.data[[level]] %in% keep_groups) %>%
    mutate(!!level := factor(.data[[level]], levels = keep_groups))

  p_base <- ggplot(
    df_level,
    aes(x = Known_Log_Lifespan_years,
        y = Predicted_Log_Lifespan_years,
        color = .data[[level]], group = .data[[level]])
  ) +
    geom_point(alpha = 0.65, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.7, color = "black") +
    coord_equal(xlim = c(X_MIN, X_MAX), ylim = c(Y_MIN, Y_MAX), expand = FALSE) +
    scale_x_continuous(breaks = BRKS, limits = c(X_MIN, X_MAX)) +
    scale_y_continuous(breaks = BRKS, limits = c(Y_MIN, Y_MAX)) +
    labs(
      x = "ln(Known lifespan years)",
      y = "ln(Predicted lifespan years)",
      color = level
    ) +
    theme_all_black +
    theme(
      axis.title.x = element_text(size = 22, color = "black"),
      axis.title.y = element_text(size = 22, color = "black"),
      axis.text.x  = element_text(size = 20, color = "black"),
      axis.text.y  = element_text(size = 20, color = "black"),
      legend.title = element_text(size = 16, color = "black"),
      legend.text  = element_text(size = 14, color = "black")
    )

  # Legend as a separate graphic
  leg_plot <- extract_legend(p_base, ncol = 2, inner_left_pt = 32, outer_left_pt = 36)

  # Save main plot without legend
  p_main <- p_base + theme(legend.position = "none")
  out_plot <- paste0("scatter_known_vs_pred_colored_", tolower(level), "_plot_ln_years.png")
  ggsave(out_plot, plot = p_main, width = out_width, height = out_height, dpi = dpi)
  message("Wrote: ", out_plot)

  # Save legend separately
  out_legend <- paste0("scatter_known_vs_pred_colored_", tolower(level), "_legend.png")
  ggsave(out_legend, plot = leg_plot, width = 6, height = 8, dpi = dpi, bg = "white")
  message("Wrote: ", out_legend)
}

# ====== Run for each taxonomic level ======
for (lev in taxonomic_levels) {
  plot_scatter_colored_by_group(lev, min_n = 3, max_groups = Inf,
                                out_width = 16, out_height = 10)
}

# ====== Whole-Protostomia scatter (ln-years) ======
plot_scatter_all_protostomia <- function(out_width = 10, out_height = 9, dpi = 300) {
  d <- filtered_data

  p_all <- ggplot(
    d,
    aes(x = Known_Log_Lifespan_years, y = Predicted_Log_Lifespan_years)
  ) +
    geom_point(alpha = 0.65, size = 1.8, color = "black") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2, color = "red") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.7, color = "black") +
    coord_equal(xlim = c(X_MIN, X_MAX), ylim = c(Y_MIN, Y_MAX), expand = FALSE) +
    scale_x_continuous(breaks = BRKS, limits = c(X_MIN, X_MAX)) +
    scale_y_continuous(breaks = BRKS, limits = c(Y_MIN, Y_MAX)) +
    labs(
      x = "ln(Known lifespan years)",
      y = "ln(Predicted lifespan years)"
    ) +
    theme_all_black +
    theme(legend.position = "none")

  # annotation (robust to fixed limits)
  fit <- lm(Predicted_Log_Lifespan_years ~ Known_Log_Lifespan_years, data = d)
  s  <- unname(coef(fit)["Known_Log_Lifespan_years"])
  r2 <- summary(fit)$r.squared
  n  <- nrow(d)
  p_all <- p_all +
    annotate("text",
             x = X_MIN + 0.02 * (X_MAX - X_MIN),
             y = Y_MAX - 0.02 * (Y_MAX - Y_MIN),
             label = sprintf("slope = %.3f, R^2 = %.3f, n = %d", s, r2, n),
             hjust = 0, vjust = 1, size = 5, color = "black")

  ggsave("scatter_known_vs_pred_all_protostomia_plot_ln_years.png",
         plot = p_all, width = out_width, height = out_height, dpi = dpi)
  message("Wrote: scatter_known_vs_pred_all_protostomia_plot_ln_years.png")
}

# Run it
plot_scatter_all_protostomia()


# ===================== 6-panel figure (ln-years, fixed axes) =====================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(cowplot)
})

# ---- Axis constants ----
X_MIN <- -4; X_MAX <- 6
Y_MIN <- -4; Y_MAX <- 6
BRKS  <- seq(-4, 6, by = 1)

# ---- Panel theme (no per-panel axis titles) ----
panel_theme <- theme_all_black +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, color = "black")
  )

# ---- Helper: one panel for a given taxonomic level ----
scatter_plot_for_level <- function(level, min_n = 3, max_groups = Inf) {
  df0 <- filtered_data %>% filter(!is.na(.data[[level]]))

  groups <- df0 %>%
    count(across(all_of(level)), name = "n_group") %>%
    arrange(desc(n_group)) %>%
    filter(n_group >= min_n)

  if (nrow(groups) == 0) return(ggplot() + theme_void())

  if (is.finite(max_groups) && nrow(groups) > max_groups) {
    groups <- groups %>% slice_head(n = max_groups)
  }
  keep <- groups %>% pull(!!sym(level))

  df <- df0 %>%
    filter(.data[[level]] %in% keep) %>%
    mutate(!!level := factor(.data[[level]], levels = keep))

  ggplot(
    df,
    aes(Known_Log_Lifespan_years, Predicted_Log_Lifespan_years,
        color = .data[[level]], group = .data[[level]])
  ) +
    geom_point(alpha = 0.65, size = 1.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.7, color = "black") +
    coord_equal(xlim = c(X_MIN, X_MAX), ylim = c(Y_MIN, Y_MAX), expand = FALSE) +
    scale_x_continuous(breaks = BRKS, limits = c(X_MIN, X_MAX)) +
    scale_y_continuous(breaks = BRKS, limits = c(Y_MIN, Y_MAX)) +
    panel_theme
}

# ---- Protostomia panel (all species pooled) ----
scatter_plot_protostomia <- function() {
  d <- filtered_data
  ggplot(d, aes(Known_Log_Lifespan_years, Predicted_Log_Lifespan_years)) +
    geom_point(alpha = 0.65, size = 1.6, color = "black") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.2, color = "red") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.7, color = "black") +
    coord_equal(xlim = c(X_MIN, X_MAX), ylim = c(Y_MIN, Y_MAX), expand = FALSE) +
    scale_x_continuous(breaks = BRKS, limits = c(X_MIN, X_MAX)) +
    scale_y_continuous(breaks = BRKS, limits = c(Y_MIN, Y_MAX)) +
    panel_theme
}

# ---- Build the six panels ----
p_proto  <- scatter_plot_protostomia() + ggtitle("Protostomia")
p_phylum <- scatter_plot_for_level("Phylum")    + ggtitle("Phylum")
p_class  <- scatter_plot_for_level("Class")     + ggtitle("Class")
p_order  <- scatter_plot_for_level("Order")     + ggtitle("Order")
p_family <- scatter_plot_for_level("Family")    + ggtitle("Family")
p_genus  <- scatter_plot_for_level("Genus")     + ggtitle("Genus")

title_theme <- theme(
  plot.title = element_text(size = 18, face = "plain", color = "black", hjust = 0.5)
)
p_proto  <- p_proto  + title_theme
p_phylum <- p_phylum + title_theme
p_class  <- p_class  + title_theme
p_order  <- p_order  + title_theme
p_family <- p_family + title_theme
p_genus  <- p_genus  + title_theme

# ---- Arrange 2×3 with a spacer column ----
left_col  <- (p_proto / p_class / p_family)
right_col <- (p_phylum / p_order / p_genus)

combined <- (left_col | plot_spacer() | right_col) +
  plot_layout(widths = c(1, 0.10, 1))


# ---- Global axis labels in dedicated white gutters & save ----
# Reserve left (for y label) and bottom (for x label) gutters by shrinking the panel area.
# All coordinates here are in [0,1] relative to the whole figure.
final <- cowplot::ggdraw() +
  # place the 2×3 panel grid, leaving room: left=7%, bottom=9%
  cowplot::draw_plot(combined, x = 0.07, y = 0.09, width = 0.93, height = 0.91) +

  # x-axis label in bottom gutter
  cowplot::draw_label(
    "ln(Known lifespan years)",
    x = 0.53, y = 0.035,        # safely inside the bottom gutter
    vjust = 0, size = 20
  ) +

  # y-axis label in left gutter
  cowplot::draw_label(
    "ln(Predicted lifespan years)",
    x = 0.02, y = 0.54,         # safely inside the left gutter
    angle = 90, vjust = 1, size = 20
  ) +

  # ensure the background is white (including gutters)
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(
  "scatter_known_vs_pred_6panel_2col3row_common_axes_ln_years.png",
  plot = final, width = 14, height = 16, dpi = 300, bg = "white"
)
