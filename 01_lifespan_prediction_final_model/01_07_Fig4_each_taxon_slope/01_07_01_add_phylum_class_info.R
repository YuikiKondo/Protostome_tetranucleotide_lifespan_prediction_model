setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_final_model/01_07_Fig4_each_taxon_slope")

# Load necessary library
library(dplyr)

# Step 1: Load both CSV files
tetranuc_data <- read.csv("filtered_species_all_lifecycle_categories_each_region_wide_tetranuc_rates_Github_check.csv", header = TRUE)
lifespan_predictions <- read.csv("log_lifespan_predictions.csv", header = TRUE)

# Step 2: Select relevant columns for merging
taxonomy_info <- tetranuc_data %>%
  select(`Organism.Name`, Phylum, Class, Order, Family, Genus, Parthenogenesis, `Fission.Budding`, `Highly.dormant.eggs`, `Dormant.state.in.other.than.egg.stage`, `More.than.one.type.of.life.cycles`, Eusociality, `Eusociality.level`, `Host.number`)

# Step 3: Merge based on matching organism name
merged_data <- lifespan_predictions %>%
  left_join(taxonomy_info, by = c("Organism_Name" = "Organism.Name"))

# Step 4: Reorder columns to put Phylum and Class first
final_data <- merged_data %>%
  select(Phylum, Class, Order, Family, Genus, everything())

# Step 5: Write to a new CSV file
write.csv(final_data, "log_lifespan_predictions_with_taxonomy.csv", row.names = FALSE)







# Load necessary library
library(dplyr)

# Step 1: Read the prediction results with taxonomy
predictions <- read.csv("log_lifespan_predictions_with_taxonomy.csv")

# Step 2: Add back-calculated lifespan in years
predictions <- predictions %>%
  mutate(
    Known_Lifespan_years = exp(Known_Log_Lifespan) / 365,
    Predicted_Lifespan_years = exp(Predicted_Log_Lifespan_Model_Averaging) / 365
  )

# Step 3: Get all dataset groups and phyla
all_groups <- unique(predictions$Dataset)
all_phyla <- unique(predictions$Phylum)

# Step 4: Define function to compute error metrics
compute_error_metrics <- function(group_name, phylum_name, df) {
  group_df <- df %>%
    filter(Dataset == group_name, if (phylum_name != "ALL") Phylum == phylum_name else TRUE)
  
  # Skip if no species
  if (nrow(group_df) == 0) return(NULL)
  
  # Log-scale metrics
  mse_log <- mean((group_df$Known_Log_Lifespan - group_df$Predicted_Log_Lifespan_Model_Averaging)^2, na.rm = TRUE)
  mae_log <- mean(abs(group_df$Known_Log_Lifespan - group_df$Predicted_Log_Lifespan_Model_Averaging), na.rm = TRUE)
  mean_err_log <- mean(group_df$Predicted_Log_Lifespan_Model_Averaging - group_df$Known_Log_Lifespan, na.rm = TRUE)
  
  # Years-scale metrics
  mse_years <- mean((group_df$Known_Lifespan_years - group_df$Predicted_Lifespan_years)^2, na.rm = TRUE)
  mae_years <- mean(abs(group_df$Known_Lifespan_years - group_df$Predicted_Lifespan_years), na.rm = TRUE)
  mean_err_years <- mean(group_df$Predicted_Lifespan_years - group_df$Known_Lifespan_years, na.rm = TRUE)
  max_err_years <- max(abs(group_df$Predicted_Lifespan_years - group_df$Known_Lifespan_years), na.rm = TRUE)
  
  # Over/Underestimated species count
  n_over <- sum(group_df$Predicted_Lifespan_years > group_df$Known_Lifespan_years, na.rm = TRUE)
  n_under <- sum(group_df$Predicted_Lifespan_years < group_df$Known_Lifespan_years, na.rm = TRUE)
  n_species <- nrow(group_df)
  proportion_over <- ifelse(n_species > 0, n_over / n_species, NA)
  
  # Create result row
  data.frame(
    Dataset = group_name,
    Phylum = phylum_name,
    Num_Species = n_species,
    N_Overestimated = n_over,
    N_Underestimated = n_under,
    Proportion_Overestimated = round(proportion_over, 3),
    MSE_Ln_Lifespan_DAYS = mse_log,
    MAE_Ln_Lifespan_DAYS = mae_log,
    Mean_Error_Ln_Lifespan_DAYS = mean_err_log,
    MSE_Lifespan_Years = mse_years,
    MAE_Lifespan_Years = mae_years,
    Mean_Error_Lifespan_Years = mean_err_years,
    Max_Error_Lifespan_Years = max_err_years
  )
}

# Step 5: Apply to all Dataset × Phylum and Dataset × ALL
error_metrics_list <- list()

for (group in all_groups) {
  # For each Phylum
  for (phylum in all_phyla) {
    result <- compute_error_metrics(group, phylum, predictions)
    if (!is.null(result)) {
      error_metrics_list[[length(error_metrics_list) + 1]] <- result
    }
  }
  # Also compute overall metrics (Phylum = "ALL")
  overall_result <- compute_error_metrics(group, "ALL", predictions)
  if (!is.null(overall_result)) {
    error_metrics_list[[length(error_metrics_list) + 1]] <- overall_result
  }
}

# Step 6: Combine and save
error_metrics_combined <- do.call(rbind, error_metrics_list)
write.csv(error_metrics_combined, "lifespan_prediction_error_metrics_by_phylum_and_all.csv", row.names = FALSE)
