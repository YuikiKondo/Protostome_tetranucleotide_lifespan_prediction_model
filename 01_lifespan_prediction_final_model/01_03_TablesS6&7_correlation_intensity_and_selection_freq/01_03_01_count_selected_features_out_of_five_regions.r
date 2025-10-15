# ====== Set working directory ======
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_and_error/01_03_TablesS6&7_correlation_intensity_and_selection_freq")

# ====== Libraries ======
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
})

# ====== Parameters ======
THRESH <- 0.0  # "selected" if Selection_Frequency > THRESH
REGIONS <- c("genome", "upstream", "exon", "intron", "downstream")

# Generate the full 256 tetranucleotides
acgt <- c("A","C","G","T")
ALL_TETRAS <- apply(expand.grid(acgt, acgt, acgt, acgt), 1, paste0, collapse = "")

# ====== Load ======
infile <- "combined_model_weights.csv"
df_raw <- read_csv(infile, show_col_types = FALSE)

# ====== Parse features & normalize regions robustly ======
normalize_region <- function(x) {
  x_up <- toupper(x)
  case_when(
    str_detect(x_up, "UPSTREAM")                      ~ "upstream",
    str_detect(x_up, "DOWNSTREAM")                    ~ "downstream",
    str_detect(x_up, "^INTRON(S)?$")                  ~ "intron",
    str_detect(x_up, "^EXON(S)?$")                    ~ "exon",
    str_detect(x_up, "^(GENOME|ALL|WHOLE)$")          ~ "genome",
    # Sometimes suffixes like EXONS_1ST or INTRON_3 may appear; strip trailing junk:
    str_detect(x_up, "^EXON")                         ~ "exon",
    str_detect(x_up, "^INTRON")                       ~ "intron",
    TRUE ~ NA_character_
  )
}

df <- df_raw %>%
  filter(Feature_Name != "(Intercept)") %>%
  mutate(
    Tetranuc = str_extract(Feature_Name, "^[ACGT]{4}"),          # ensure exactly 4-mer
    Region_raw = str_replace(Feature_Name, "^[ACGT]+_Rate_", ""),
    Region = normalize_region(Region_raw),
    Selection_Frequency = coalesce(Selection_Frequency, 0)       # treat NA as 0
  ) %>%
  filter(!is.na(Tetranuc), !is.na(Region), Region %in% REGIONS)

# ====== Flag selection per (Tetranuc, Region) with NA-safe logic ======
df_sel <- df %>%
  mutate(selected = Selection_Frequency > THRESH) %>%
  group_by(Tetranuc, Region) %>%
  summarise(selected = any(selected, na.rm = TRUE), .groups = "drop")

# ====== Fill to full 256 × 5 grid (missing pairs -> selected = FALSE) ======
df_sel_complete <- df_sel %>%
  complete(Tetranuc = ALL_TETRAS, Region = REGIONS, fill = list(selected = FALSE))

# ====== Count selected regions per tetranucleotide ======
per_tetranuc_counts <- df_sel_complete %>%
  group_by(Tetranuc) %>%
  summarise(n_regions_selected = sum(selected, na.rm = TRUE), .groups = "drop")

# ====== Summarize distribution across 0..5 with correct proportions ======
summary_base <- per_tetranuc_counts %>%
  count(n_regions_selected, name = "n_tetranucs") %>%
  complete(n_regions_selected = 0:5, fill = list(n_tetranucs = 0L)) %>%
  arrange(n_regions_selected)

total_n <- sum(summary_base$n_tetranucs)

summary_counts <- summary_base %>%
  mutate(
    label = case_when(
      n_regions_selected == 0 ~ "none of 5 regions",
      n_regions_selected == 1 ~ "only 1 region",
      n_regions_selected == 2 ~ "2 regions",
      n_regions_selected == 3 ~ "3 regions",
      n_regions_selected == 4 ~ "4 regions",
      n_regions_selected == 5 ~ "all 5 regions"
    ),
    prop = if (total_n > 0) n_tetranucs / total_n else NA_real_
  )

# ====== Save ======
out_prefix <- if (THRESH == 0) "gt0" else paste0("gt", THRESH)
write_csv(per_tetranuc_counts, paste0("tetranuc_n_regions_selected_", out_prefix, ".csv"))
write_csv(summary_counts,     paste0("tetranuc_n_regions_selected_summary_", out_prefix, ".csv"))

# ====== Console sanity checks ======
cat("Unique tetranucs (should be 256): ", n_distinct(per_tetranuc_counts$Tetranuc), "\n")
bad_regions <- df_raw %>%
  transmute(Region_raw = str_replace(Feature_Name, "^[ACGT]+_Rate_", "")) %>%
  mutate(Norm = normalize_region(Region_raw)) %>%
  filter(!is.na(Region_raw) & is.na(Norm)) %>%
  distinct(Region_raw)
if (nrow(bad_regions) > 0) {
  cat("Unmapped region strings (inspect to extend normalize_region()):\n")
  print(bad_regions, n = Inf)
}
print(summary_counts, n = Inf)
cat("Check sum(prop) =", sum(summary_counts$prop, na.rm = TRUE), "\n")
