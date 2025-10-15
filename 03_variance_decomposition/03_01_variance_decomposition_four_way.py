# save as: compute_r2_and_partition_4way.py
import pandas as pd
import numpy as np
from itertools import chain, combinations
from pathlib import Path

# ----------------------- CONFIG: filenames -----------------------
# Keys are subset labels of {T, M, P, G}
FILES = {
    "T":   "log_lifespan_predictions_tetra.csv",
    "M":   "log_lifespan_predictions_mechanics.csv",
    "P":   "log_lifespan_predictions_phylo.csv",
    "G":   "log_lifespan_predictions_genomesize.csv",

    "TM":  "log_lifespan_predictions_tetra_mechanics.csv",
    "TP":  "log_lifespan_predictions_tetra_phylo.csv",
    "TG":  "log_lifespan_predictions_tetra_genomesize.csv",
    "MP":  "log_lifespan_predictions_mechanics_phylo.csv",
    "MG":  "log_lifespan_predictions_mechanics_genomesize.csv",
    "PG":  "log_lifespan_predictions_phylo_genomesize.csv",

    "TMP": "log_lifespan_predictions_tetra_mechanics_phylo.csv",
    "TMG": "log_lifespan_predictions_tetra_mechanics_genomesize.csv",
    "TPG": "log_lifespan_predictions_tetra_phylo_genomesize.csv",
    "MPG": "log_lifespan_predictions_mechanics_phylo_genomesize.csv",

    "TMPG":"log_lifespan_predictions_tetra_mechanics_phylo_genomesize.csv",
}

SETS = ["T", "M", "P", "G"]  # order we’ll use consistently

# ----------------------- helpers -----------------------
def powerset(iterable):
    """all non-empty subsets, as tuple of labels, sorted by label for consistency"""
    s = list(iterable)
    return [tuple(x) for r in range(1, len(s)+1) for x in combinations(s, r)]

def label_of(subset):
    """subset ('M','T') -> canonical label 'MT' sorted by SETS order"""
    return "".join(sorted(subset, key=SETS.index))

def r2_from_file(path: Path) -> float:
    """
    Compute R^2 on the TEST subset using:
        R^2 = 1 - sum((y - yhat)^2) / sum((y - y_mean)^2)
    Requires columns: Dataset, Known_Log_Lifespan, Predicted_Log_Lifespan_Model_Averaging
    """
    df = pd.read_csv(path)
    if "Dataset" not in df.columns:
        raise ValueError(f"{path} is missing a 'Dataset' column.")
    test = df[df["Dataset"].astype(str).str.lower() == "test"].copy()
    if test.empty:
        raise ValueError(f"{path} has no Test rows.")
    y = test["Known_Log_Lifespan"].to_numpy()
    yhat = test["Predicted_Log_Lifespan_Model_Averaging"].to_numpy()
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot

# ----------------------- 1) Compute R² for all 15 models -----------------------
all_models = [label_of(s) for s in powerset(SETS)]   # 15 labels in canonical order
missing = [lab for lab in all_models if lab not in FILES]
if missing:
    raise ValueError(f"FILES dict is missing entries for: {missing}")

r2 = {}
for lab in all_models:
    p = Path(FILES[lab])
    if not p.exists():
        raise FileNotFoundError(f"Missing file: {p.resolve()}")
    r2[lab] = r2_from_file(p)

r2_table = pd.DataFrame({"Model": all_models, "R2": [r2[m] for m in all_models]})
r2_table.to_csv("R2_all_models_4way.csv", index=False)

# ----------------------- 2) Build & solve the 4-way partition -----------------------
# Atoms are also all non-empty subsets of {T,M,P,G} (15 atoms).
atoms = [label_of(s) for s in powerset(SETS)]  # same list as models, but now interpreted as "regions"

# For a model with predictor set S, it can explain any atom J that includes
# at least one predictor from S (i.e., J ∩ S ≠ ∅). That’s the standard
# “attribution to a set” definition used in ecological variation partitioning.
def subset_from_label(lab):  # 'TMP' -> ('T','M','P')
    return tuple(l for l in SETS if l in set(lab))

def can_explain(model_lab, atom_lab):
    S = set(subset_from_label(model_lab))
    J = set(subset_from_label(atom_lab))
    return int(len(S.intersection(J)) > 0)

# Build linear system A x = b
A = np.zeros((len(all_models), len(atoms)), dtype=float)
for i, mlab in enumerate(all_models):
    for j, alab in enumerate(atoms):
        A[i, j] = can_explain(mlab, alab)

b = np.array([r2[m] for m in all_models], dtype=float)

# Solve for atom values (unique & all overlaps)
# (Matrix is full rank for 4 sets; numerical noise can create very small negatives)
x = np.linalg.solve(A, b)

partition_df = pd.DataFrame({
    "Atom": atoms,                 # e.g., 'T', 'TM', 'TMP', 'TMPG', etc.
    "R2_amount": x
}).sort_values("Atom", key=lambda s: s.str.len())  # shorter labels first
total_full = r2["TMPG"]
partition_df["Percent_of_R2_TMPG"] = 100 * partition_df["R2_amount"] / total_full
partition_df.to_csv("variance_partitioning_T_M_P_G_atoms.csv", index=False)

# ----------------------- 3) Exact “inside T” decomposition (8 parts) -----------------------
# Inside T are the atoms that include 'T':
inside_T_labels = [a for a in atoms if "T" in a]
inside_T_values = {lab: x[atoms.index(lab)] for lab in inside_T_labels}
R2_T = r2["T"]

# Break out into the 8 disjoint parts that sum to R2_T:
# T, TM, TP, TG, TMP, TMG, TPG, TMPG
inside_T_exact = pd.DataFrame({
    "Within_T_component_exact": [
        "Unique_T_only", "Overlap_T∩M_only", "Overlap_T∩P_only", "Overlap_T∩G_only",
        "Overlap_T∩M∩P", "Overlap_T∩M∩G", "Overlap_T∩P∩G", "Overlap_T∩M∩P∩G"
    ],
    "Atom_Label": ["T", "TM", "TP", "TG", "TMP", "TMG", "TPG", "TMPG"],
    "R2_amount": [
        inside_T_values.get("T", 0.0),
        inside_T_values.get("TM", 0.0),
        inside_T_values.get("TP", 0.0),
        inside_T_values.get("TG", 0.0),
        inside_T_values.get("TMP", 0.0),
        inside_T_values.get("TMG", 0.0),
        inside_T_values.get("TPG", 0.0),
        inside_T_values.get("TMPG", 0.0),
    ],
})
inside_T_exact["Percent_of_R2_T"] = 100 * inside_T_exact["R2_amount"] / R2_T
inside_T_exact.to_csv("inside_T_partition_exact_4way.csv", index=False)

# ----------------------- 4) Simple order-specific unique ΔR² (added last) -----------------------
# “What does each set add when added last to the other three?”
unique_last = pd.DataFrame({
    "Set": ["T","M","P","G"],
    "Delta_R2_added_last": [
        r2["TMPG"] - r2["MPG"],
        r2["TMPG"] - r2["TPG"],
        r2["TMPG"] - r2["TMG"],
        r2["TMPG"] - r2["TMP"],
    ]
})
unique_last.to_csv("unique_last_deltas.csv", index=False)

# ----------------------- 5) Sanity checks (printed) -----------------------
print("Saved:")
print(" - R2_all_models_4way.csv")
print(" - variance_partitioning_T_M_P_G_atoms.csv")
print(" - inside_T_partition_exact_4way.csv")
print(" - unique_last_deltas.csv")
print("\nSanity checks:")
print(f"Sum of all 15 atoms = {partition_df['R2_amount'].sum():.6f} (should equal R2_TMPG = {total_full:.6f})")
print(f"Sum inside T (8 parts) = {inside_T_exact['R2_amount'].sum():.6f} (should equal R2_T = {R2_T:.6f})")
