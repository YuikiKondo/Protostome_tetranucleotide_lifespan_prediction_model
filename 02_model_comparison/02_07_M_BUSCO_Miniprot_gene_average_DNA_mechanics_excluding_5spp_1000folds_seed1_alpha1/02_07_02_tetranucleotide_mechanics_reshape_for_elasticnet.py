import pandas as pd

# Load the CSV file
df = pd.read_csv('/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/di_tri_tetra_mechanics_excluding_5spp_Github_check/filtered_species_with_all_indexes_di_tri_tetra_combined_Github_check.csv')

# Define the columns to retain and their order
columns_to_keep = [
    'Assembly Identifier', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Organism Name', 'Average_lifespan_days', 'Region',
    "DNA_bendability", "X_Displacement", "Y_Displacement", "Inclination", "Tip", "Shear", "Stretch", "Stagger", "Buckle", "Propel", "Opening", "Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Chi", "Phase", "Amplitude", "Hydrogen_bond", "Stacking_energy", "Solvation", "Shift", "Slide", "Rise", "Tilt", "Roll", "Twist", "A_philicity_energy"
]

# Select and rearrange the columns
df_filtered = df[columns_to_keep]

# Save the filtered data to a new CSV file
df_filtered.to_csv('/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/di_tri_tetra_mechanics_excluding_5spp_Github_check/elasticnet_all_genes_average/filtered_species_with_zero_lifecycle_6_categories_rearranged_di_tri_tetra_mechanics.csv', index=False)

print("Filtered and rearranged data saved as 'filtered_species_with_zero_lifecycle_6_categories_rearranged_di_tri_tetra_mechanics.csv'")
