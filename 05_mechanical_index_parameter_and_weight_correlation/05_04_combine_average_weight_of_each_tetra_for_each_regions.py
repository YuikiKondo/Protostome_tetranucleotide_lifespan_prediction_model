import pandas as pd
import re

# ======================
# File paths (adjust if needed)
# ======================
mechanics_file = "tetranucleotide_mechanics_full_zscore.csv"
weights_file = "combined_model_weights.csv"
output_file = "tetranucleotide_mechanics_full_zscore_with_weights_by_region.csv"

# ======================
# Load input files
# ======================
mechanics_df = pd.read_csv(mechanics_file)
weights_df = pd.read_csv(weights_file)

# ======================
# Parse tetranucleotide and region from Feature_Name
# ======================
weights_df = weights_df.copy()
weights_df["Tetranuc"] = weights_df["Feature_Name"].str.extract(r"^([ACGT]{4})", expand=False)

# Grab the last underscore-delimited token as region
weights_df["Region"] = weights_df["Feature_Name"].str.extract(
    r"(upstream\d+|downstream\d+|exon|intron|genome)$", expand=False
)

# ✅ Rename upstream2→upstream, downstream2→downstream
weights_df["Region"] = weights_df["Region"].replace({
    "upstream2": "upstream",
    "downstream2": "downstream"
})

# Keep only rows where both parsed fields exist
weights_df = weights_df.dropna(subset=["Tetranuc", "Region"])

# ======================
# 1) Average across *all* regions
# ======================
avg_across = (
    weights_df.groupby("Tetranuc", as_index=False)["Average_Weight"]
    .mean()
    .rename(columns={"Average_Weight": "Average_Weight_Tetranuc"})
)

# ======================
# 2) Average *per region* (one column per region)
# ======================
region_avg_wide = (
    weights_df.groupby(["Tetranuc", "Region"])["Average_Weight"]
    .mean()
    .unstack("Region")
)

# Prefix columns for clarity
region_avg_wide = region_avg_wide.add_prefix("Average_Weight_").reset_index()

# Ensure a consistent column order
desired_regions = ["exon", "intron", "upstream", "downstream", "genome"]
ordered_cols = ["Tetranuc"] + [
    f"Average_Weight_{r}" for r in desired_regions if f"Average_Weight_{r}" in region_avg_wide.columns
] + [c for c in region_avg_wide.columns if c not in ["Tetranuc"] + [f"Average_Weight_{r}" for r in desired_regions]]
region_avg_wide = region_avg_wide[ordered_cols]

# ======================
# 3) Combine across- and per-region averages
# ======================
weights_summary = avg_across.merge(region_avg_wide, on="Tetranuc", how="left")

# ======================
# 4) Merge with mechanics dataframe on motif
# ======================
merged_df = mechanics_df.merge(weights_summary, left_on="Motif", right_on="Tetranuc", how="left")
merged_df = merged_df.drop(columns=["Tetranuc"])

# ======================
# 5) Save output
# ======================
merged_df.to_csv(output_file, index=False)
print(f"✅ Done! Saved file with per-region and across-region weights: {output_file}")
