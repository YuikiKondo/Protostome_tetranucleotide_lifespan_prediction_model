import pandas as pd

# ======================
# File paths (adjust if needed)
# ======================
mechanics_file = "tetranucleotide_mechanics_full_zscore.csv"
weights_file = "combined_model_weights.csv"
output_file = "tetranucleotide_mechanics_full_zscore_with_weights.csv"

# ======================
# Load input files
# ======================
mechanics_df = pd.read_csv(mechanics_file)
weights_df = pd.read_csv(weights_file)

# ======================
# Extract tetranucleotide from Feature_Name
# ======================
# Example: "AAAA_Rate_downstream2" → "AAAA"
weights_df["Tetranuc"] = weights_df["Feature_Name"].str.split("_").str[0]

# ======================
# Compute average weight per tetranucleotide
# ======================
avg_weights = (
    weights_df.groupby("Tetranuc")["Average_Weight"]
    .mean()
    .reset_index()
    .rename(columns={"Average_Weight": "Average_Weight_Tetranuc"})
)

# ======================
# Merge with mechanics dataframe
# ======================
# mechanics_df["Motif"] contains tetranucleotides
merged_df = mechanics_df.merge(avg_weights, left_on="Motif", right_on="Tetranuc", how="left")

# Drop the extra key column
merged_df = merged_df.drop(columns=["Tetranuc"])

# ======================
# Save output
# ======================
merged_df.to_csv(output_file, index=False)

print(f"✅ Done! Saved file with weights: {output_file}")
