import pandas as pd

# Load the CSV file
df = pd.read_csv('/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/dinucleotide_count_excluding_5spp/filtered_species_with_zero_lifecycle_6_categories_Github_check.csv')

# Define the columns to retain and their order
columns_to_keep = [
    'Assembly Identifier', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Organism Name', 'Average_lifespan_days', 'Region',
    'AA_Rate', 'AT_Rate', 'AG_Rate', 'AC_Rate', 'TA_Rate', 'TT_Rate', 'TG_Rate', 'TC_Rate',
    'GA_Rate', 'GT_Rate', 'GG_Rate', 'GC_Rate', 'CA_Rate', 'CT_Rate', 'CG_Rate', 'CC_Rate'
]

# Select and rearrange the columns
df_filtered = df[columns_to_keep]

# Save the filtered data to a new CSV file
df_filtered.to_csv('/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/dinucleotide_count_excluding_5spp/elasticnet_all_gene_average/filtered_species_with_zero_lifecycle_6_categories_rearranged_dinuc_rates_Github_check.csv', index=False)

print("Filtered and rearranged data saved as 'filtered_species_with_zero_lifecycle_6_categories_rearranged_dinuc_rates_Github_check.csv'")
