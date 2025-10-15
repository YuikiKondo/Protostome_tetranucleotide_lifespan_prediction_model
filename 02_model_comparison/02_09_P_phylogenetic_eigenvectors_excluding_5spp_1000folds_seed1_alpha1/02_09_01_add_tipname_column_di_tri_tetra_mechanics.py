import pandas as pd

# Define file paths 
# Add tip names corresponding to the phylogenetic tree that will be used for PGLS
busco_file = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_exon_intron_separated/gene_region_separated_sequences_no_internal_stop_3divisible/tetranucleotide_count/DNA_mechanics/PGLS/busco_species_tree_tip_alternative_names.csv"

# filtered_species_file
filtered_species_file = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/di_tri_tetra_mechanics_excluding_5spp_Github_check/filtered_species_with_all_indexes_di_tri_tetra_combined_Github_check.csv"

output_file = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/di_tri_tetra_mechanics_excluding_5spp_Github_check/elasticnet_all_genes_average_phylogenetic_matrix/filtered_species_with_tip_name_di_tri_tetra_mechanics_Github_check.csv"

# Load the CSV files
busco_df = pd.read_csv(busco_file)
filtered_species_df = pd.read_csv(filtered_species_file)

# Merge based on "Organism Name" in filtered_species_df and "Species" in busco_df
merged_df = filtered_species_df.merge(busco_df[['Species', 'Tip_Name']], 
                                      left_on="Organism Name", 
                                      right_on="Species", 
                                      how="left")

# Drop the redundant "Species" column (optional)
merged_df.drop(columns=["Species"], inplace=True)

# Save the new file
merged_df.to_csv(output_file, index=False)

print(f"File saved successfully: {output_file}")
