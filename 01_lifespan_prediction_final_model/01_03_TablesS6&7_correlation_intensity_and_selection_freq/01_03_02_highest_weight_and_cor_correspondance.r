# ====== Set working directory ======
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_and_error/01_03_TablesS6&7_correlation_intensity_and_selection_freq")

# ====== Libraries ======
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

# ---------- Input ----------
infile <- "combined_model_weights.csv"
df_raw <- read_csv(infile, show_col_types = FALSE)

# ---------- Parse + prep ----------
df2 <- df_raw %>%
  filter(Feature_Name != "(Intercept)") %>%
  mutate(
    Tetranuc = str_extract(Feature_Name, "^[ACGT]{4}"),
    Region   = str_replace(Feature_Name, "^[ACGT]+_Rate_", ""),
    Selection_Frequency = coalesce(Selection_Frequency, 0),
    abs_weight = abs(Average_Weight),
    abs_corr   = abs(Lifespan_Feature_Correlation_All_Species)
  ) %>%
  filter(!is.na(Tetranuc), !is.na(Region))

# ---------- Helper: top region(s) by a metric (handles ties) ----------
top_by_metric <- function(data, metric_col, metric_name){
  data %>%
    group_by(Tetranuc) %>%
    filter(.data[[metric_col]] == max(.data[[metric_col]], na.rm = TRUE)) %>%
    summarise(
      !!paste0("max_", metric_name) := first(.data[[metric_col]]),
      !!paste0("regions_top_", metric_name) := paste(sort(unique(Region)), collapse = "|"),
      !!paste0("n_top_", metric_name, "_ties") := n(),
      .groups = "drop"
    )
}

# ---------- One run for a given strict threshold (> thr) ----------
run_for_threshold <- function(thr){
  d <- df2 %>% filter(Selection_Frequency > thr)
  total_features <- nrow(d)  # <-- NEW: total number of features passing the threshold

  # If no features pass, return empty but with total_features = 0
  if (total_features == 0) {
    comp <- tibble(
      Tetranuc = character(),
      max_abs_weight = numeric(),
      regions_top_abs_weight = character(),
      n_top_abs_weight_ties = integer(),
      max_abs_corr = numeric(),
      regions_top_abs_corr = character(),
      n_top_abs_corr_ties = integer(),
      any_match = logical(),
      strict_equal = logical()
    )
    summary_counts <- tibble(
      threshold = thr,
      total_features = 0L,         # <-- NEW COLUMN
      n_tetranuc = 0L,
      n_any_match = 0L,
      n_strict_equal = 0L,
      prop_any_match = NA_real_,
      prop_strict_equal = NA_real_
    )
    return(list(comp = comp, summary = summary_counts))
  }

  top_w <- top_by_metric(d, "abs_weight", "abs_weight")
  top_c <- top_by_metric(d, "abs_corr",  "abs_corr")

  comp <- top_w %>%
    inner_join(top_c, by = "Tetranuc") %>%
    mutate(
      any_match = map2_lgl(strsplit(regions_top_abs_weight, "\\|"),
                           strsplit(regions_top_abs_corr,  "\\|"),
                           ~ length(intersect(.x, .y)) > 0),
      strict_equal = (n_top_abs_weight_ties == 1 &
                      n_top_abs_corr_ties  == 1 &
                      regions_top_abs_weight == regions_top_abs_corr)
    )

  summary_counts <- comp %>%
    summarise(
      n_tetranuc        = n(),
      n_any_match       = sum(any_match),
      n_strict_equal    = sum(strict_equal),
      prop_any_match    = ifelse(n_tetranuc > 0, n_any_match / n_tetranuc, NA_real_),
      prop_strict_equal = ifelse(n_tetranuc > 0, n_strict_equal / n_tetranuc, NA_real_)
    ) %>%
    mutate(
      threshold = thr,
      total_features = total_features,  # <-- NEW COLUMN
      .before = 1
    )

  list(comp = comp, summary = summary_counts)
}

# ---------- Thresholds to evaluate (> thr) ----------
thresholds <- c(0.50, 0.30, 0.10, 0.00)

# Nice labels for filenames
thr_label <- function(thr) paste0("gt", format(round(100*thr, 0), trim = TRUE), "pct")

# ---------- Run and write files ----------
all_summaries <- map_dfr(thresholds, function(thr){
  res <- run_for_threshold(thr)

  # per-threshold outputs
  write_csv(res$comp,    paste0("tetranuc_weight_vs_corr_selFreq_", thr_label(thr), ".csv"))
  write_csv(res$summary, paste0("tetranuc_weight_vs_corr_summary_selFreq_", thr_label(thr), ".csv"))

  res$summary
})

# ---------- Combined summary (now includes total_features) ----------
write_csv(all_summaries, "tetranuc_weight_vs_corr_summary_ALL_thresholds.csv")
print(all_summaries, n = Inf)
