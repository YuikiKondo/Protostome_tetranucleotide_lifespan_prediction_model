import pandas as pd

# Load the CSV file
df = pd.read_csv("/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/trinucleotide_count/elasticnet_all_gene_average_Github_check/filtered_species_with_zero_lifecycle_6_categories_rearranged_augustus_all_gene_trinuc_rates_Github_check.csv")

# Fill fill 0 in empty (NaN) cells in the target column
df['Average_lifespan_days'] = df['Average_lifespan_days'].fillna(0)

# Identify relevant columns
key_columns = ['Assembly Identifier', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Organism Name', 'Average_lifespan_days', 'Region']
trinuc_columns = [col for col in df.columns if col not in key_columns]

# Reshape the dataframe
reshaped_df = df.melt(id_vars=key_columns, value_vars=trinuc_columns, var_name='Trinucleotide', value_name='Rate')
reshaped_df['Trinucleotide'] = reshaped_df['Trinucleotide'] + "_" + reshaped_df['Region']
reshaped_df = reshaped_df.drop(columns=['Region'])

# Pivot the table so that each unique Organism Name has all trinucleotide rates as separate columns
final_df = reshaped_df.pivot_table(index=['Assembly Identifier', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Organism Name', 'Average_lifespan_days'], columns='Trinucleotide', values='Rate').reset_index()

# Save the reshaped dataframe
final_df.to_csv("filtered_species_with_zero_lifecycle_6_categories_each_region_wide_augustus_all_gene_trinuc_rates_Github_check.csv", index=False)

print("Reshaped file saved as filtered_species_with_zero_lifecycle_6_categories_each_region_wide_augustus_all_gene_trinuc_rates_Github_check.csv")
