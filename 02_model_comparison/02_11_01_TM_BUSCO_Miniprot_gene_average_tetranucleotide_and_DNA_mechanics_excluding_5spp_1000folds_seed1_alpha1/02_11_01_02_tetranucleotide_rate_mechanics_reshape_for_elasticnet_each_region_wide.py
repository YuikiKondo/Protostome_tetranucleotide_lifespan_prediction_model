import pandas as pd

# Load the CSV file
df = pd.read_csv("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/di_tri_tetra_mechanics_excluding_5spp_Github_check/elasticnet_all_genes_average_tetra_rate_and_30_mechanics/filtered_species_with_zero_lifecycle_6_categories_rearranged_di_tri_tetra_mechanics.csv")

# Identify relevant columns
key_columns = ['Assembly Identifier', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Organism Name', 'Average_lifespan_days', 'Region']
rate_mechanics_columns = [col for col in df.columns if col not in key_columns]

# Reshape the dataframe
reshaped_df = df.melt(id_vars=key_columns, value_vars=rate_mechanics_columns, var_name='Rate_Mechanics', value_name='Index')
reshaped_df['Rate_Mechanics'] = reshaped_df['Rate_Mechanics'] + "_" + reshaped_df['Region']
reshaped_df = reshaped_df.drop(columns=['Region'])

# Pivot the table so that each unique Organism Name has all mechanical indices as separate columns
final_df = reshaped_df.pivot_table(index=['Assembly Identifier', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Organism Name', 'Average_lifespan_days'], columns='Rate_Mechanics', values='Index').reset_index()

# Save the reshaped dataframe
final_df.to_csv("filtered_species_with_zero_lifecycle_6_categories_each_region_wide_tetra_rate_di_tri_tetra_mechanics.csv", index=False)

print("Reshaped file saved as filtered_species_with_zero_lifecycle_6_categories_each_region_wide_tetra_rate_di_tri_tetra_mechanics.csv")
