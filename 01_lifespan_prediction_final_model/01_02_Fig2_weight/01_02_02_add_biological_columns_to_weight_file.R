setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_and_error/01_02_Fig2_weight")

library(dplyr)
library(stringr)
library(readr)

# Load data
df <- read_csv("combined_model_weights.csv")

# Remove intercept row if present
df <- df %>% filter(Feature_Name != "(Intercept)")

# Extract tetranucleotide (letters at the start) and region (after "_Rate_")
df <- df %>%
  mutate(
    Feature_Name = toupper(Feature_Name),
    Tetranuc = str_extract(Feature_Name, "^[ACGT]+"),          # restrict to canonical bases
    Region_raw = str_extract(Feature_Name, "(?<=_RATE_).*$")   # grab region suffix (now uppercased)
  )

# Normalize Region to canonical buckets
df <- df %>%
  mutate(
    Region_norm = case_when(
      is.na(Region_raw)                        ~ NA_character_,
      str_detect(Region_raw, "UPSTREAM")       ~ "upstream",
      str_detect(Region_raw, "DOWNSTREAM")     ~ "downstream",
      str_detect(Region_raw, "^INTRON$")       ~ "intron",
      str_detect(Region_raw, "^EXON$")         ~ "exon",
      str_detect(Region_raw, "^GENOME$")       ~ "genome",
      TRUE                                     ~ tolower(Region_raw)
    )
  )

# Precompute base counts and length
df <- df %>%
  mutate(
    n_len = nchar(Tetranuc),
    n_GC  = str_count(Tetranuc, "[GC]"),
    n_AT  = str_count(Tetranuc, "[AT]"),
    n_pur = str_count(Tetranuc, "[AG]"),
    n_pyr = str_count(Tetranuc, "[CT]")
  )

# Add requested motif flags
df <- df %>%
  mutate(
    CpG_including    = as.integer(str_detect(Tetranuc, "CG")),
    CAA_including    = as.integer(str_detect(Tetranuc, "CAA")),   # NEW column
    GC_rich          = as.integer(n_len > 0 & (n_GC / n_len) >= 0.75),
    AT_rich          = as.integer(n_len > 0 & (n_AT / n_len) >= 0.75),
    Purine_rich      = as.integer(n_len > 0 & (n_pur / n_len) >= 0.75),
    Pyrimidine_rich  = as.integer(n_len > 0 & (n_pyr / n_len) >= 0.75),
    ATGG_related     = as.integer(str_detect(Tetranuc, "ATG|TGG")),
    Unclassified = as.integer(
      (CpG_including + CAA_including + GC_rich + AT_rich +
         Purine_rich + Pyrimidine_rich + ATGG_related) == 0
    )
  )

# Region one-hot columns
df <- df %>%
  mutate(
    Exon       = as.integer(Region_norm == "exon"),
    Intron     = as.integer(Region_norm == "intron"),
    Upstream   = as.integer(str_detect(coalesce(Region_norm, ""), "^upstream")),
    Downstream = as.integer(str_detect(coalesce(Region_norm, ""), "^downstream")),
    Genome     = as.integer(Region_norm == "genome")
  )

# Tidy up columns: keep both raw and normalized region if useful
df_out <- df %>%
  select(
    everything(),
    -n_len, -n_GC, -n_AT, -n_pur, -n_pyr
  ) %>%
  rename(Region = Region_norm)

# Save new file
write_csv(df_out, "combined_model_weights_with_bio.csv")
message("Saved file: combined_model_weights_with_bio.csv")
