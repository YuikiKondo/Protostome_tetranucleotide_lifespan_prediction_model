setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_and_error/01_06_TableS9_clean_tetra_columns_Miniprot")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ===== Inputs / Outputs =====
train_in  <- "data_definitive_train.csv"
test_in   <- "data_definitive_test.csv"

train_out <- "data_definitive_train_filtered.csv"
test_out  <- "data_definitive_test_filtered.csv"
combo_out <- "data_definitive_combined_filtered.csv"

# Columns to remove (keep Organism.Name)
cols_unnecessary <- c(
  "Phylum", "Class", "Order", "Family", "Genus",
  "Lifespan", "log_lifespan"
)

# ---- helper: drop unnecessary columns, remove *1 columns, rename *2 -> base ----
process_file_remove <- function(path) {
  df <- read_csv(path, show_col_types = FALSE) %>%
    select(-any_of(cols_unnecessary)) %>%              # drop listed metadata cols (keep Organism.Name)
    select(-matches("downstream1|upstream1"))          # drop *1 columns

  # rename *2 -> base
  nm <- names(df)
  nm <- gsub("downstream2", "downstream", nm, fixed = TRUE)
  nm <- gsub("upstream2",   "upstream",   nm, fixed = TRUE)
  names(df) <- nm

  df
}

# ---- train ----
df_train <- process_file_remove(train_in)
write_csv(df_train, train_out)

# ---- test ----
df_test <- process_file_remove(test_in)
write_csv(df_test, test_out)

# ---- combine with consistent column order ----
id_cols   <- c("Assembly Identifier", "Organism.Name")
all_cols  <- union(names(df_train), names(df_test))
rest_cols <- setdiff(all_cols, id_cols)
order_cols <- c(id_cols, sort(rest_cols))

df_train <- df_train %>% select(any_of(order_cols))
df_test  <- df_test  %>% select(any_of(order_cols))

df_combined <- bind_rows(df_train, df_test)
write_csv(df_combined, combo_out)

message("Wrote files:\n- ", train_out, "\n- ", test_out, "\n- ", combo_out)
