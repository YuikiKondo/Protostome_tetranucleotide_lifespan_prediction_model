import pandas as pd
from sklearn.preprocessing import StandardScaler

# === Input/Output files ===
input_file = "tetranucleotide_mechanics_full.csv"
output_file = "tetranucleotide_mechanics_full_zscore.csv"

# === Load ===
df = pd.read_csv(input_file)

# === Identify feature columns (everything except 'Motif') ===
features = [c for c in df.columns if c != "Motif"]

# === Z-score transform ===
scaler = StandardScaler()
z_vals = scaler.fit_transform(df[features])

# Create new DataFrame with z-scored values
df_z = pd.DataFrame(z_vals, columns=[f"{c}_z" for c in features])
df_z.insert(0, "Motif", df["Motif"])

# === Save ===
df_z.to_csv(output_file, index=False)

print(f"Z-scored file saved: {output_file}, shape = {df_z.shape}")
