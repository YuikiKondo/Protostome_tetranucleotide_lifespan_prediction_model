setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_and_error/01_03_FigS2_frequently_selected_feature_dinuc_trinucleotides")

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
    Tetranuc     = str_extract(Feature_Name, "^[ACGT]+"),
    Region_raw   = str_extract(Feature_Name, "(?<=_RATE_).*$")
  )

# Normalize Region to canonical buckets
df <- df %>%
  mutate(
    Region_norm = case_when(
      is.na(Region_raw)                    ~ NA_character_,
      str_detect(Region_raw, "UPSTREAM")   ~ "upstream",
      str_detect(Region_raw, "DOWNSTREAM") ~ "downstream",
      str_detect(Region_raw, "^INTRON$")   ~ "intron",
      str_detect(Region_raw, "^EXON$")     ~ "exon",
      str_detect(Region_raw, "^GENOME$")   ~ "genome",
      TRUE                                 ~ tolower(Region_raw)
    )
  )

# ---- Richness flags (keep these) ----
df <- df %>%
  mutate(
    n_len = nchar(Tetranuc),
    n_GC  = str_count(Tetranuc, "[GC]"),
    n_AT  = str_count(Tetranuc, "[AT]"),
    n_pur = str_count(Tetranuc, "[AG]"),
    n_pyr = str_count(Tetranuc, "[CT]"),
    GC_rich         = as.integer(n_len > 0 & (n_GC / n_len) >= 0.75),
    AT_rich         = as.integer(n_len > 0 & (n_AT / n_len) >= 0.75),
    Purine_rich     = as.integer(n_len > 0 & (n_pur / n_len) >= 0.75),
    Pyrimidine_rich = as.integer(n_len > 0 & (n_pyr / n_len) >= 0.75),
    AT_exclusive    = as.integer(n_len > 0 & (n_AT / n_len) == 1)
  )

# ---- 16 dinucleotide flags ----
dinucs <- c("AA","AC","AG","AT","CA","CC","CG","CT",
            "GA","GC","GG","GT","TA","TC","TG","TT")

dinuc_cols <- setNames(
  lapply(dinucs, function(d) as.integer(str_detect(df$Tetranuc, d))),
  paste0(dinucs, "_including")
)

df <- bind_cols(df, as.data.frame(dinuc_cols))

# ---- 64 trinucleotide flags ----
bases   <- c("A","C","G","T")
trinucs <- as.vector(outer(outer(bases, bases, paste0), bases, paste0))  # AAA..TTT

trinuc_cols <- setNames(
  lapply(trinucs, function(tr) as.integer(str_detect(df$Tetranuc, tr))),
  paste0(trinucs, "_including")
)

df <- bind_cols(df, as.data.frame(trinuc_cols))

# ---- Region one-hot columns ----
df <- df %>%
  mutate(
    Exon       = as.integer(Region_norm == "exon"),
    Intron     = as.integer(Region_norm == "intron"),
    Upstream   = as.integer(str_detect(coalesce(Region_norm, ""), "^upstream")),
    Downstream = as.integer(str_detect(coalesce(Region_norm, ""), "^downstream")),
    Genome     = as.integer(Region_norm == "genome")
  )

# ---- Output (drop intermediate counts) ----
df_out <- df %>%
  select(-n_len, -n_GC, -n_AT, -n_pur, -n_pyr) %>%
  rename(Region = Region_norm)

write_csv(df_out, "combined_model_weights_with_bio.csv")
message("Saved file: combined_model_weights_with_bio.csv (with richness + 16 dinucs + 64 trinucs)")
