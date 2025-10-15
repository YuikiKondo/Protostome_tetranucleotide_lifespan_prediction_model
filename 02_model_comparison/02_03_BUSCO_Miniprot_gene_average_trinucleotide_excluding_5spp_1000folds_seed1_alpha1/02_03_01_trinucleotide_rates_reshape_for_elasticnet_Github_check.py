import pandas as pd

# Load the CSV file
df = pd.read_csv('/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/trinucleotide_count_excluding_5spp/filtered_species_with_zero_lifecycle_6_categories_Github_check.csv')

# Define the columns to retain and their order
columns_to_keep = [
    'Assembly Identifier', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Organism Name', 'Average_lifespan_days', 'Region', 
    "AAA_Rate", "AAT_Rate", "AAG_Rate", "AAC_Rate", "ATA_Rate", "ATT_Rate", "ATG_Rate", "ATC_Rate", "AGA_Rate", "AGT_Rate", "AGG_Rate", "AGC_Rate", "ACA_Rate", "ACT_Rate", "ACG_Rate", "ACC_Rate", "TAA_Rate", "TAT_Rate", "TAG_Rate", "TAC_Rate", "TTA_Rate", "TTT_Rate", "TTG_Rate", "TTC_Rate", "TGA_Rate", "TGT_Rate", "TGG_Rate", "TGC_Rate", "TCA_Rate", "TCT_Rate", "TCG_Rate", "TCC_Rate", "GAA_Rate", "GAT_Rate", "GAG_Rate", "GAC_Rate", "GTA_Rate", "GTT_Rate", "GTG_Rate", "GTC_Rate", "GGA_Rate", "GGT_Rate", "GGG_Rate", "GGC_Rate", "GCA_Rate", "GCT_Rate", "GCG_Rate", "GCC_Rate", "CAA_Rate", "CAT_Rate", "CAG_Rate", "CAC_Rate", "CTA_Rate", "CTT_Rate", "CTG_Rate", "CTC_Rate", "CGA_Rate", "CGT_Rate", "CGG_Rate", "CGC_Rate", "CCA_Rate", "CCT_Rate", "CCG_Rate", "CCC_Rate"
]

# Select and rearrange the columns
df_filtered = df[columns_to_keep]

# Save the filtered data to a new CSV file
df_filtered.to_csv('/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/trinucleotide_count_excluding_5spp/elasticnet_all_gene_average/filtered_species_with_zero_lifecycle_6_categories_rearranged_trinuc_rates_Github_check.csv', index=False)

print("Filtered and rearranged data saved as 'filtered_species_with_zero_lifecycle_6_categories_rearranged_trinuc_rates_Github_check.csv'")
