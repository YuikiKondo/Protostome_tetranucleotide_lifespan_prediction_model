setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/06_cross_clade_models_by_Ecdysozoa_and_Spiralia/06_03_Trian_by_Spiralia_Test_by_Ecdysozoa")

# Load libraries
library(dplyr)
library(readr)

# Load data
data <- read_csv("filtered_species_with_zero_lifecycle_6_categories_each_region_wide_tetranuc_rates_Github_check.csv")

# Keep only rows where exon columns are not empty
data <- data %>%
  filter(if_all(matches("exon"), ~ . != ""))

# Rename columns
data <- data %>%
  rename(
    Organism.Name = `Organism Name`,
    Lifespan = `Average_lifespan_days`
  )

# Order data alphabetically
data <- data[order(data$Organism.Name), ]

# Log-transform lifespan
data <- data %>%
  mutate(
    log_lifespan = ifelse(Lifespan > 0, log(Lifespan), NA)
  )

# Assign phyla to Spiralia or Ecdysozoa
spiralia_phyla <- c("Mollusca", "Annelida", "Platyhelminthes", "Brachiopoda", "Bryozoa")
ecdysozoa_phyla <- c("Nematoda", "Arthropoda (Chelicerata)", "Arthropoda (Mandibulata)", "Acanthocephala")

# Split train/test by phylum
train_data <- data %>% filter(Phylum %in% spiralia_phyla)
test_data  <- data %>% filter(Phylum %in% ecdysozoa_phyla)

# Save outputs
write_csv(train_data, "data_definitive_train.csv")
write_csv(test_data, "data_definitive_test.csv")

# Quick summary
cat("Total species:", nrow(data), "\n")
cat("Train (Spiralia):", nrow(train_data), "\n")
cat("Test (Ecdysozoa):", nrow(test_data), "\n")
cat("Train phyla:", unique(train_data$Phylum), "\n")
cat("Test phyla:", unique(test_data$Phylum), "\n")
