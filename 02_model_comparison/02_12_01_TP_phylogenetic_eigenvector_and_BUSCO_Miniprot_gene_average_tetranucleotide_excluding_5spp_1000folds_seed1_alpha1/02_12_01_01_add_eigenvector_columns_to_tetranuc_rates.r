# Set working directory if needed
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/02_model_comparison/02_12_01_TP_phylogenetic_eigenvector_and_BUSCO_Miniprot_gene_average_tetranucleotide_excluding_5spp_1000folds_seed1_alpha1")

# Load required libraries
library(dplyr)
library(readr)

# Load your main tetranucleotide file
tetranuc_df <- read_csv("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/03_BUSCO_Miniprot_gene_average_tetranucleotide_excluding_5spp_100folds_seed1/filtered_species_with_zero_lifecycle_6_categories_each_region_wide_tetranuc_rates_Github_check.csv")

# Load eigenvector file
eigen_df <- read_csv("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/05_phylogenetic_correlation_matrix_excluding_5spp_100folds_seed1/phylo_eigenvectors_with_lifespan.csv")

# Keep just the join key + eigenvectors from the eigen file
eigen_df_filtered <- eigen_df %>%
  select(Assembly.Identifier, starts_with("eigenvector"))

# Merge by Assembly ID (space vs dot)
merged_df <- tetranuc_df %>%
  left_join(eigen_df_filtered, by = c("Assembly Identifier" = "Assembly.Identifier"))


# (Optional) quick checks
cat("Rows in tetranucleotide rates:", nrow(tetranuc_df), "\n")
cat("Rows in merged   :", nrow(merged_df), "\n")
cat("New columns added:", sum(grepl("^eigenvector_", names(merged_df))), "\n")


# Save as new CSV
write_csv(merged_df, "filtered_species_with_zero_lifecycle_6_categories_each_region_wide_tetranuc_rates_with_eigenvectors_Github_check.csv")
