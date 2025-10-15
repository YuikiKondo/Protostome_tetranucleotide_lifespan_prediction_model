# Set working directory if needed
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/02_14_01_MP_phylogenetic_eigenvector_and_BUSCO_Miniprot_gene_average_DNA_mechanics_excluding_5spp_1000folds_seed1_alpha1")

# Load required libraries
library(dplyr)
library(readr)

# Load your main mechanics file
mechanics_df <- read_csv("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/04_01_M_BUSCO_Miniprot_gene_average_DNA_mechanics_excluding_5spp_1000folds_seed1_alpha1/filtered_species_with_zero_lifecycle_6_categories_each_region_wide_di_tri_tetra_mechanics.csv")

# Load eigenvector file
eigen_df <- read_csv("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/05_phylogenetic_correlation_matrix_excluding_5spp_1000folds_seed1_alpha1/phylo_eigenvectors_with_lifespan.csv")

# Keep just the join key + eigenvectors from the eigen file
eigen_df_filtered <- eigen_df %>%
  select(Assembly.Identifier, starts_with("eigenvector"))

# Merge by Assembly ID (space vs dot)
merged_df <- mechanics_df %>%
  left_join(eigen_df_filtered, by = c("Assembly Identifier" = "Assembly.Identifier"))


# (Optional) quick checks
cat("Rows in mechanics:", nrow(mechanics_df), "\n")
cat("Rows in merged   :", nrow(merged_df), "\n")
cat("New columns added:", sum(grepl("^eigenvector_", names(merged_df))), "\n")


# Save as new CSV
write_csv(merged_df, "filtered_species_with_zero_lifecycle_6_categories_each_region_wide_di_tri_tetra_mechanics_with_eigenvectors.csv")
