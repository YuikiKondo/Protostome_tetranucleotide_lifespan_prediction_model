import pandas as pd
import itertools

# ==== File paths (change to your own) ====
di_file = "Sharma_Supplemantary_file_dinucleotide.csv"
tri_file = "Sharma_Supplemantary_file_trinuc_revised_20250919.csv"
tetra_file = "Sharma_Supplemantary_file_tetranucleotide.csv"
output_file = "tetranucleotide_mechanics_full.csv"

# ==== Load data ====
di_df = pd.read_csv(di_file).rename(columns={"Dinucleotide": "Motif"})
tri_df = pd.read_csv(tri_file).rename(columns={"Trinucleotide": "Motif"})
tetra_df = pd.read_csv(tetra_file).rename(columns={"Tetranucleotide": "Motif"})

# ==== Generate all 256 tetranucleotides ====
bases = ["A", "C", "G", "T"]
tetranucleotides = ["".join(p) for p in itertools.product(bases, repeat=4)]

# ==== Helper functions ====
def get_dinucleotides(tetra):
    return [tetra[0:2], tetra[1:3], tetra[2:4]]

def get_trinucleotides(tetra):
    return [tetra[0:3], tetra[1:4]]

# ==== Lift dinucleotide indices to tetranucleotide level ====
di_features = [c for c in di_df.columns if c != "Motif"]
di_dict = di_df.set_index("Motif").to_dict(orient="index")

tetra_di = []
for tet in tetranucleotides:
    sub_di = get_dinucleotides(tet)
    values = {}
    for f in di_features:
        values[f] = sum(di_dict[sub][f] for sub in sub_di) / len(sub_di)
    tetra_di.append({"Motif": tet, **values})
tetra_di_df = pd.DataFrame(tetra_di)

# ==== Lift trinucleotide indices to tetranucleotide level ====
tri_features = [c for c in tri_df.columns if c != "Motif"]
tri_dict = tri_df.set_index("Motif").to_dict(orient="index")

tetra_tri = []
for tet in tetranucleotides:
    sub_tri = get_trinucleotides(tet)
    values = {}
    for f in tri_features:
        values[f] = sum(tri_dict[sub][f] for sub in sub_tri) / len(sub_tri)
    tetra_tri.append({"Motif": tet, **values})
tetra_tri_df = pd.DataFrame(tetra_tri)

# ==== Merge all three sources ====
merged = tetra_df.merge(tetra_di_df, on="Motif").merge(tetra_tri_df, on="Motif")

# ==== Save to CSV ====
merged.to_csv(output_file, index=False)

print(f"Saved combined mechanics table with shape {merged.shape} to {output_file}")
