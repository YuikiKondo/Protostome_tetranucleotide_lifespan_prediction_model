# Set working directory if needed
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/02_model_comparison/02_21_02_TPMG_genome_size_phylogenetic_eigenvector_and_Augustus_all_gene_average_tetranucleotide_and_DNA_mechanics_excluding_5spp_1000folds_seed1_alpha1")

# Load libraries
library(dplyr)
library(readr)
library(caret)

# Step 1: Load the data
data <- read_csv("filtered_species_with_zero_lifecycle_6_categories_each_region_wide_tetra_rate_and_di_tri_tetra_mechanics_with_eigenvectors_and_genome_size.csv")

# Step 2: Rename columns for convenience
data <- data %>%
  rename(
    Organism.Name = `Organism Name`,
    Lifespan = `Average_lifespan_days`
  )

# Step 2.1: Order data alphabetically by Organism.Name
data <- data[order(data$Organism.Name), ]

# Step 3: Log-transform lifespan
data <- data %>%
  mutate(
    log_lifespan = ifelse(Lifespan > 0, log(Lifespan), NA)
  )

# Step 4: Filter definitive life cycle species (input file already includes only definitive species)
data_definitive <- data 

# Step 5: 70/30 split within each Order
set.seed(1)
train_list <- list()
test_list <- list()

for (ord in unique(data_definitive$Order)) {
  subset_data <- data_definitive %>% filter(Order == ord)
  n <- nrow(subset_data)
  
  if (n == 2) {
    idx <- sample(1:2, 1)
    train_list[[ord]] <- subset_data[idx, , drop = FALSE]
    test_list[[ord]] <- subset_data[-idx, , drop = FALSE]
  } else if (n >= 3) {
    desired_n_train <- round(n * 0.7)
    # use p = 0.9 here (a high proportion) to give createDataPartition() enough flexibility to form balanced strata.
    # this may return more rows than desired_n_train, which leads to the next step.
    idx <- createDataPartition(subset_data$log_lifespan, times = 1, p = 0.9, list = FALSE)
    
    # Downsample to exact number if needed
    # If too many samples were selected by createDataPartition(), downsample to exactly desired_n_train using sample().
    if (length(idx) > desired_n_train) {
      idx <- sample(idx, desired_n_train)
    }

    train_list[[ord]] <- subset_data[idx, ]
    test_list[[ord]] <- subset_data[-idx, ]
  }
}

# Identify singleton orders
singleton_orders <- data_definitive %>%
  group_by(Order) %>%
  summarise(n_species = n()) %>%
  filter(n_species == 1) %>%
  pull(Order)

# Extract all singleton-species rows
singleton_species <- data_definitive %>% filter(Order %in% singleton_orders)

# Split 70/30 at species level (not per order since there's only 1 per order)
set.seed(1)
train_singleton_idx <- createDataPartition(singleton_species$log_lifespan, p = 0.7, list = FALSE)

train_singleton <- singleton_species[train_singleton_idx, ]
test_singleton <- singleton_species[-train_singleton_idx, ]

# Append singleton splits to the existing lists
train_data <- bind_rows(bind_rows(train_list), train_singleton)
test_data <- bind_rows(bind_rows(test_list), test_singleton)


# Step 6: Save output
write_csv(train_data, "data_definitive_train.csv")
write_csv(test_data, "data_definitive_test.csv")

nrow(data_definitive)
nrow(train_data)
nrow(test_data)


train_data %>% filter(Order == "Diptera") %>% nrow()
test_data %>% filter(Order == "Diptera") %>% nrow()

train_data %>% filter(Order == "Monhysterida") %>% nrow()
test_data %>% filter(Order == "Monhysterida") %>% nrow()

train_data %>% filter(Order == "Trombidiformes") %>% nrow()
test_data %>% filter(Order == "Trombidiformes") %>% nrow()
