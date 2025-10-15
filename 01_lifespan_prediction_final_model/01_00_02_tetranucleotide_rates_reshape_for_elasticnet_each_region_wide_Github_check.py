import pandas as pd

# Load the CSV file
df = pd.read_csv("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/tetranucleotide_count/elasticnet_all_gene_average_6categories_test/filtered_species_with_zero_lifecycle_6_categories_rearranged_tetranuc_rates_Github_check.csv")

# Fill fill 0 in empty (NaN) cells in the target column
df['Average_lifespan_days'] = df['Average_lifespan_days'].fillna(0)

# Identify relevant columns
key_columns = ['Assembly Identifier', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Organism Name', 'Average_lifespan_days', 'Region']
tetranuc_columns = [col for col in df.columns if col not in key_columns]

# Reshape the dataframe
reshaped_df = df.melt(id_vars=key_columns, value_vars=tetranuc_columns, var_name='Tetranucleotide', value_name='Rate')
reshaped_df['Tetranucleotide'] = reshaped_df['Tetranucleotide'] + "_" + reshaped_df['Region']
reshaped_df = reshaped_df.drop(columns=['Region'])

# Pivot the table so that each unique Organism Name has all tetranucleotide rates as separate columns
final_df = reshaped_df.pivot_table(index=['Assembly Identifier', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Organism Name', 'Average_lifespan_days'], columns='Tetranucleotide', values='Rate').reset_index()

# Save the reshaped dataframe
final_df.to_csv("filtered_species_with_zero_lifecycle_6_categories_each_region_wide_tetranuc_rates_Github_check.csv", index=False)

print("Reshaped file saved as filtered_species_with_zero_lifecycle_6_categories_each_region_wide_tetranuc_rates_Github_check.csv")
