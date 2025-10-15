import pandas as pd

# Load the three CSV files
tetranuc_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/tetranucleotide_count/DNA_mechanics/filtered_species_with_all_indexes_Github_check.csv"
trinuc_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/trinucleotide_count/DNA_mechanics/filtered_species_with_all_indexes_Github_check.csv"
dinuc_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/dinucleotide_count/DNA_mechanics/filtered_species_with_all_indexes_Github_check.csv"

df_tetranuc = pd.read_csv(tetranuc_file)
df_trinuc = pd.read_csv(trinuc_file)
df_dinuc = pd.read_csv(dinuc_file)

# Specify the columns to merge from tetranucleotide
columns_tri = ["DNA_bendability", "X_Displacement", "Y_Displacement", "Inclination", "Tip", "Shear", "Stretch", "Stagger", "Buckle", "Propel", "Opening", "Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Chi", "Phase", "Amplitude", "Hydrogen_bond", "Stacking_energy", "Solvation"]

# Merge tetranucleotide file with trinucleotide mechanics
df_merged = df_tetranuc.merge(
    df_trinuc[["Organism Name", "Region"] + columns_tri],
    on=["Organism Name", "Region"],
    how="left"
)

# Merge the dinucleotide "A_philicity_energy" column
df_merged = df_merged.merge(
    df_dinuc[["Organism Name", "Region", "A_philicity_energy"]],
    on=["Organism Name", "Region"],
    how="left"
)

# Save the merged DataFrame
output_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/di_tri_tetra_mechanics_Github_check/filtered_species_with_all_indexes_di_tri_tetra_combined_Github_check.csv"
df_merged.to_csv(output_file, index=False)
print(f"Merged file saved to: {output_file}")
