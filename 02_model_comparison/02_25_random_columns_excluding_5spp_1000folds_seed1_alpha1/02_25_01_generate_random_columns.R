setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/02_model_comparison/02_25_random_columns_excluding_5spp_1000folds_seed1_alpha1")

## ==== Load data ====
data <- read.csv("filtered_species_with_zero_lifecycle_6_categories_each_region_wide_tetranuc_rates_Github_check.csv",
                 stringsAsFactors = FALSE)

## ==== Remove columns with "_Rate" in the name ====
data_no_rate <- data[ , !grepl("_Rate", names(data))]

## ==== Generate 1280 random columns (values from 0 to 1) ====
set.seed(123)  # for reproducibility
random_matrix <- matrix(runif(nrow(data_no_rate) * 1280, min = 0, max = 1),
                        nrow = nrow(data_no_rate), ncol = 1280)

## Assign column names: Random_1 ... Random_1280
colnames(random_matrix) <- paste0("Random_", seq_len(1280))

## Convert to data.frame and combine with original
random_df <- as.data.frame(random_matrix)
final_data <- cbind(data_no_rate, random_df)

## ==== Save result ====
write.csv(final_data,
          "filtered_species_with_zero_lifecycle_6_categories_each_region_wide_random_1280columns.csv",
          row.names = FALSE)
