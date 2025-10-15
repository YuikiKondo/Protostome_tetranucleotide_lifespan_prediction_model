# ============================
# Count Proportions by Trait × Region (Selection_Frequency > 0.1)
# ============================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ---- Load data ----
df <- read_csv("combined_model_weights_with_bio.csv", show_col_types = FALSE)

traits <- c("CpG_including", "CAA_including", "ATGG_related",
            "GC_rich", "AT_rich",
            "Purine_rich", "Pyrimidine_rich", "Unclassified")

# Ensure tetranucleotide parsed
if (!"Tetranuc" %in% names(df)) {
  if (!"Feature_Name" %in% names(df)) stop("Feature_Name column not found. Needed to extract Tetranuc.")
  df <- df |>
    mutate(
      Feature_Name = toupper(Feature_Name),
      Tetranuc = str_extract(Feature_Name, "^[ACGT]+")
    )
}

# Canonical region list
all_regions <- df |> distinct(Region) |> pull(Region) |> sort()

# ---- 1) Region totals (ALL features, filtered) ----
region_totals <- df |>
  filter(Average_Weight != 0, Selection_Frequency > 0.1) |>
  group_by(Region) |>
  summarize(
    Positive_Total_Region = sum(Average_Weight > 0, na.rm = TRUE),
    Negative_Total_Region = sum(Average_Weight < 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  complete(Region = all_regions,
           fill = list(Positive_Total_Region = 0,
                       Negative_Total_Region = 0))

# ---- 2) Trait × Region counts (filtered) ----
summarize_trait_by_region <- function(trait_name) {
  df |>
    filter(.data[[trait_name]] == 1,
           Average_Weight != 0,
           Selection_Frequency > 0.1) |>
    group_by(Region) |>
    summarize(
      Positive_Count = sum(Average_Weight > 0, na.rm = TRUE),
      Negative_Count = sum(Average_Weight < 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    complete(Region = all_regions,
             fill = list(Positive_Count = 0, Negative_Count = 0)) |>
    mutate(Trait = trait_name) |>
    relocate(Trait, Region)
}

trait_region_counts <- bind_rows(lapply(traits, summarize_trait_by_region))

# ---- 3) Possible tetranucleotide counts (trait × region, ignoring filters) ----
possible_trait_region <- lapply(traits, function(tr) {
  df |>
    filter(.data[[tr]] == 1) |>
    group_by(Region) |>
    summarize(Possible_Tetranucleotide_Number = n_distinct(Tetranuc), .groups = "drop") |>
    complete(Region = all_regions, fill = list(Possible_Tetranucleotide_Number = 0)) |>
    mutate(Trait = tr) |>
    relocate(Trait, Region)
}) |> bind_rows()

# ---- 4) Compute proportions + add possible counts and expected proportion ----
trait_region_props <- trait_region_counts |>
  left_join(region_totals, by = "Region") |>
  left_join(possible_trait_region, by = c("Trait", "Region")) |>
  mutate(
    Positive_Prop = ifelse(Positive_Total_Region > 0,
                           Positive_Count / Positive_Total_Region, NA_real_),
    Negative_Prop = ifelse(Negative_Total_Region > 0,
                           Negative_Count / Negative_Total_Region, NA_real_),
    Expected_Prop = Possible_Tetranucleotide_Number / 256
  ) |>
  arrange(Trait, Region)

# ---- 5) Save outputs ----
write_csv(trait_region_props, "count_proportions_by_trait_and_region_sel0.1.csv")
write_csv(region_totals, "count_region_totals_ALL_sel0.1.csv")

message("Saved: count_proportions_by_trait_and_region_sel0.1.csv (with Possible_Tetranucleotide_Number and Expected_Prop)")
message("Saved: count_region_totals_ALL_sel0.1.csv")
