# Prepares amino acid composition data for the AA composition baseline clustering analyses.
# Reads raw DE-STRESS data, filters to the non-redundant AF2 set (pLDDT >= 70),
# extracts the 20 standard composition_* columns (composition_UNK excluded),
# scales with standard / robust / minmax, and writes scaled CSVs plus a labels file
# to data/processed_data_aa_composition/af2/{scaler}/.

# 0. Importing packages---------------------------------------------------------
import os
import sys
import pickle

import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler, RobustScaler, StandardScaler

# 1. Defining variables---------------------------------------------------------

scaling_method_list = ["standard", "robust", "minmax"]

raw_data_path = "data/raw_data/"
processed_data_path = "data/processed_data/"
output_base_path = "data/processed_data_aa_composition/af2/"

raw_destress_data_af2_path = raw_data_path + "destress_data_af2.csv"
af2_plddt_scores_path = raw_data_path + "af2_plddt_scores.csv"
af2_structures_nonredundant_path = raw_data_path + "af2_structures_non_redundant.csv"
af2_labels_nonredundant_path = processed_data_path + "af2/labels_nonredundant.csv"

# 20 standard amino acids only - composition_UNK excluded
AA_COLS = [
    "composition_ALA",
    "composition_CYS",
    "composition_ASP",
    "composition_GLU",
    "composition_PHE",
    "composition_GLY",
    "composition_HIS",
    "composition_ILE",
    "composition_LYS",
    "composition_LEU",
    "composition_MET",
    "composition_ASN",
    "composition_PRO",
    "composition_GLN",
    "composition_ARG",
    "composition_SER",
    "composition_THR",
    "composition_VAL",
    "composition_TRP",
    "composition_TYR",
]

LABEL_COLS = [
    "design_name",
    "organism_scientific_name",
    "organism_group",
    "organism_group2",
    "subcellular_location",
]

# 2. Reading and filtering raw data---------------------------------------------

print("Reading raw AF2 DE-STRESS data...")
raw_destress_data_af2 = pd.read_csv(
    raw_destress_data_af2_path,
    usecols=["design_name"] + AA_COLS,
)

# Attach mean pLDDT and filter to pLDDT >= 70
print("Attaching pLDDT scores and filtering pLDDT >= 70...")
af2_plddt_scores = pd.read_csv(af2_plddt_scores_path)
af2_plddt_scores["design_name"] = (
    af2_plddt_scores["Filename"].str.replace(".pdb", "").astype(str)
)

raw_destress_data_af2 = raw_destress_data_af2.merge(
    af2_plddt_scores[["design_name", "Mean_PLDDT"]],
    on="design_name",
    how="left",
)

raw_destress_data_af2 = raw_destress_data_af2[
    raw_destress_data_af2["Mean_PLDDT"] >= 70
].reset_index(drop=True)

# Filter to non-redundant set
print("Filtering to non-redundant set...")
af2_structures_nonredundant = pd.read_csv(af2_structures_nonredundant_path)

raw_destress_data_af2 = raw_destress_data_af2[
    raw_destress_data_af2["design_name"].isin(
        af2_structures_nonredundant["design_name"].to_list()
    )
].reset_index(drop=True)

# Join organism labels
print("Joining organism labels...")
af2_labels = pd.read_csv(af2_labels_nonredundant_path)

data_with_labels = raw_destress_data_af2.merge(
    af2_labels[LABEL_COLS],
    on="design_name",
    how="inner",
)

print(f"Rows after filtering and label join: {len(data_with_labels)}")

# Drop any rows with missing composition values
data_with_labels = data_with_labels.dropna(subset=AA_COLS).reset_index(drop=True)
print(f"Rows after dropping missing composition values: {len(data_with_labels)}")

# 3. Saving labels file---------------------------------------------------------

os.makedirs(output_base_path, exist_ok=True)

labels_out = data_with_labels[LABEL_COLS]
labels_out.to_csv(output_base_path + "labels_nonredundant.csv", index=False)
print(f"Saved labels to {output_base_path}labels_nonredundant.csv")

# 4. Scaling and saving composition features------------------------------------

composition_data = data_with_labels[AA_COLS].reset_index(drop=True)

scaler_map = {
    "standard": StandardScaler(),
    "robust": RobustScaler(),
    "minmax": MinMaxScaler(),
}

for scaling_method in scaling_method_list:
    print(f"Scaling with {scaling_method}...")

    scaler = scaler_map[scaling_method]
    scaled_array = scaler.fit_transform(composition_data)
    scaled_df = pd.DataFrame(scaled_array, columns=AA_COLS)

    output_path_scaled = output_base_path + scaling_method + "/"
    os.makedirs(output_path_scaled, exist_ok=True)

    scaled_df.to_csv(
        output_path_scaled + "aa_composition_scaled_nonredundant.csv",
        index=False,
    )

    with open(output_path_scaled + scaling_method + "_scaler.pkl", "wb") as f:
        pickle.dump(scaler, f)

    print(f"  Saved scaled data and scaler to {output_path_scaled}")

print("Done.")
