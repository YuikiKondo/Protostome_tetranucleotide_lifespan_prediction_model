# Purpose: Calculate correlation matrix and eigenvectors based on phylogenetic closeness among species & create a figure of lifespan with phylogenetic tree (optional).

setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/02_model_comparison/02_09_P_phylogenetic_eigenvectors_excluding_5spp_1000folds_seed1_alpha1")

# Load necessary libraries
library(ape)
library(phytools)

# Load the combined species data (for lifespan, organism name and Tip_Name) and the phylogenetic tree
combined_data <- read.csv("filtered_species_with_tip_name_di_tri_tetra_mechanics_Github_check.csv") #output file of 05_01_add_tipname_column_di_tri_tetra_mechanics.py.

tree <- read.tree("/Users/yuikikondo/Desktop/Protostomes_PGLS/opentree15.1_tree/labelled_supertree/labelled_supertree_ottnames.tre")

# Calculate branch lengths on the full tree
tree <- compute.brlen(tree)

# Filter for species with lifespan data
species_with_lifespan <- combined_data[!is.na(combined_data$Average_lifespan_days), ]

cat("Unique Organism.Name count:", length(unique(species_with_lifespan $Organism.Name)), "\n")
cat("Unique Tip_Name count:", length(unique(species_with_lifespan $Tip_Name)), "\n")


# Get the list of species' tip names (Tip_Name) for matching with the tree
species_tip_names <- species_with_lifespan$Tip_Name

# Prune the tree to include only the species with lifespan data
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, species_tip_names))

# Check how many species remain in the pruned tree
cat("Number of species in the pruned tree:", length(pruned_tree$tip.label), "\n")

# Identify species that were not included in the pruned tree (but were in species_tip_names)
missing_species <- setdiff(species_tip_names, pruned_tree$tip.label)
cat("Missing species:\n")
print(missing_species)

# If all species were pruned correctly, proceed with the analysis
if (length(missing_species) == 0) {
    # Calculate the correlation matrix using vcv with corr = TRUE
    correlation_matrix <- vcv(pruned_tree, corr = TRUE)
    cat("Dimensions of correlation_matrix:", dim(correlation_matrix), "\n")
    
    # Perform eigenvector decomposition on the correlation matrix
    eigen_data <- eigen(correlation_matrix)
    
    # Calculate the minimum number of components to explain at least 95% of the variance
    num_components <- sum(cumsum(eigen_data$values) / sum(eigen_data$values) < 0.95) + 1
    corr_components <- eigen_data$vectors[, 1:num_components]
    cat("Dimensions of corr_components:", dim(corr_components), "\n")
    
    # Assign row and column names for clarity
    rownames(corr_components) <- rownames(correlation_matrix)
    colnames(corr_components) <- paste0("PC", 1:ncol(corr_components))
    
    # Create a data frame with species names as the first column
    eigenvector_df <- data.frame(
        phylotree_tip_name = rownames(corr_components),
        corr_components
    )
    
    # Save the eigenvector matrix to a CSV file with species names
    write.csv(eigenvector_df, "eigenvector_matrix.csv", row.names = FALSE)
    
    # Map Tip_Name to Organism.Name for bar graph display
    lifespan_data <- setNames(species_with_lifespan$Average_lifespan_days, 
                              species_with_lifespan$Tip_Name)
    organism_names <- setNames(species_with_lifespan$Organism.Name, 
                               species_with_lifespan$Tip_Name)
    
    # Optional: Plot the pruned tree with a bar graph of lifespans beside it
    par(mfrow = c(1, 2))  # Set up side-by-side layout

    # Plot the pruned tree with smaller tip labels
    png("pruned_tree_lifespan.png", width = 1200, height = 600, res = 150)  # Adjust resolution
    par(mfrow = c(1, 2))       # Set up side-by-side layout

    plot(pruned_tree, cex = 0.5)
    title("Pruned Tree with Species Having Lifespan Data")
    
    # Plot the bar graph with Organism.Name beside the tree
    barplot(lifespan_data[pruned_tree$tip.label], horiz = TRUE, las = 1,
            names.arg = organism_names[pruned_tree$tip.label],
            main = "Lifespan of Species", xlab = "Lifespan (days)")
    dev.off() 
    
    cat("Eigenvector matrix saved as eigenvector_matrix.csv\n")
    cat("Number of components retained to explain 95% of variance:", num_components, "\n")

} else {
    cat("There are species missing from the pruned tree.\n")
}
