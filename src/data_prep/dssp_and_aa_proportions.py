import pandas as pd

AA_COLS = [
    ("composition_ALA", "A"),
    ("composition_CYS", "C"),
    ("composition_ASP", "D"),
    ("composition_GLU", "E"),
    ("composition_PHE", "F"),
    ("composition_GLY", "G"),
    ("composition_HIS", "H"),
    ("composition_ILE", "I"),
    ("composition_LYS", "K"),
    ("composition_LEU", "L"),
    ("composition_MET", "M"),
    ("composition_ASN", "N"),
    ("composition_PRO", "P"),
    ("composition_GLN", "Q"),
    ("composition_ARG", "R"),
    ("composition_SER", "S"),
    ("composition_THR", "T"),
    ("composition_VAL", "V"),
    ("composition_TRP", "W"),
    ("composition_TYR", "Y"),
]

SS_COLS = [
    ("ss_prop_alpha_helix",   "Alpha Helix"),
    ("ss_prop_beta_bridge",   "Beta Bridge"),
    ("ss_prop_beta_strand",   "Beta Strand"),
    ("ss_prop_3_10_helix",    "3 10 Helix"),
    ("ss_prop_pi_helix",      "Pi Helix"),
    ("ss_prop_hbonded_turn",  "Hbond Turn"),
    ("ss_prop_bend",          "Bend"),
    ("ss_prop_loop",          "Loop"),
]

dssp_output_path  = "analysis/data_exploration/dssp_bin_proportions.csv"
ss_res_output_path = "analysis/data_exploration/ss_residue_proportions.csv"
aa_output_path    = "analysis/data_exploration/aa_proportions.csv"

datasets = [
    ("AFDB pLDDT > 70",                    "data/processed_data/af2/labels.csv"),
    ("AFDB pLDDT > 70 & non-redundant",    "data/processed_data/af2/labels_nonredundant.csv"),
    ("PDB",                                "data/processed_data/pdb/labels.csv"),
]

raw_destress_paths = {
    "AFDB pLDDT > 70":                  "data/raw_data/destress_data_af2.csv",
    "AFDB pLDDT > 70 & non-redundant":  "data/raw_data/destress_data_af2.csv",
    "PDB":                              "data/raw_data/destress_data_pdb_082024.csv",
}

dssp_row_order = ["Alpha Helix", "Beta Bridge", "Beta Strand", "3 10 Helix", "Pi Helix", "Hbond Turn", "Bend", "Loop", "Mixed"]

dssp_data  = {}
ss_res_data = {}
aa_data    = {}

SS_PROP_COLS = [col for col, _ in SS_COLS]
AA_PROP_COLS = [col for col, _ in AA_COLS]

for dataset_name, path in datasets:
    df = pd.read_csv(path)

    # Join ss_prop_* columns from raw DEStress data
    raw = pd.read_csv(raw_destress_paths[dataset_name], usecols=["design_name"] + SS_PROP_COLS + AA_PROP_COLS)
    df  = df.merge(raw, on="design_name", how="left")

    # --- DSSP bin proportions (per protein) ---
    counts = df["dssp_bin"].value_counts()
    props  = df["dssp_bin"].value_counts(normalize=True).mul(100)
    for bin_label in counts.index:
        if bin_label not in dssp_data:
            dssp_data[bin_label] = {}
        dssp_data[bin_label][dataset_name] = f"{props[bin_label]:.1f}% ({counts[bin_label]:,})"

    # --- SS residue proportions (length-weighted across all sequences) ---
    lengths = df["full_sequence"].dropna().str.len()
    valid   = df.loc[lengths.index]
    for col, label in SS_COLS:
        if col in valid.columns:
            mask = valid[col].notna()
            col_lengths = lengths[mask]
            weighted_pct = (valid.loc[mask, col] * col_lengths).sum() / col_lengths.sum() * 100 if col_lengths.sum() > 0 else 0.0
        else:
            weighted_pct = 0.0
        if label not in ss_res_data:
            ss_res_data[label] = {}
        ss_res_data[label][dataset_name] = f"{weighted_pct:.1f}%"

    # --- AA identity proportions (length-weighted, standard AAs only) ---
    for col, label in AA_COLS:
        if col in valid.columns:
            mask = valid[col].notna()
            col_lengths = lengths[mask]
            weighted_pct = (valid.loc[mask, col] * col_lengths).sum() / col_lengths.sum() * 100 if col_lengths.sum() > 0 else 0.0
        else:
            weighted_pct = 0.0
        if label not in aa_data:
            aa_data[label] = {}
        aa_data[label][dataset_name] = f"{weighted_pct:.2f}%"

# --- Write DSSP bin table ---
ordered_bins = list(dssp_row_order)
ordered_bins += [b for b in dssp_data if b not in dssp_row_order]

dssp_rows = []
for bin_label in ordered_bins:
    row = {"Secondary structure (per protein)": bin_label}
    for dataset_name, _ in datasets:
        row[dataset_name] = dssp_data.get(bin_label, {}).get(dataset_name, "0.0% (0)")
    dssp_rows.append(row)

dssp_result = pd.DataFrame(dssp_rows)
dssp_result.to_csv(dssp_output_path, index=False)
print("=== DSSP bin proportions (per protein) ===")
print(dssp_result.to_string(index=False))
print(f"\nSaved to {dssp_output_path}\n")

# --- Write SS residue table ---
ss_res_rows = []
for _, label in SS_COLS:
    row = {"Secondary structure (per residue)": label}
    for dataset_name, _ in datasets:
        row[dataset_name] = ss_res_data[label].get(dataset_name, "0.0%")
    ss_res_rows.append(row)

ss_res_result = pd.DataFrame(ss_res_rows)
ss_res_result.to_csv(ss_res_output_path, index=False)
print("=== Secondary structure proportions (per residue) ===")
print(ss_res_result.to_string(index=False))
print(f"\nSaved to {ss_res_output_path}\n")

# --- Write AA table ---
aa_rows = []
for _, label in AA_COLS:
    row = {"Standard amino acid (per residue)": label}
    for dataset_name, _ in datasets:
        row[dataset_name] = aa_data[label].get(dataset_name, "0.00%")
    aa_rows.append(row)

aa_result = pd.DataFrame(aa_rows)
aa_result.to_csv(aa_output_path, index=False)
print("=== Standard amino acid proportions (per residue) ===")
print(aa_result.to_string(index=False))
print(f"\nSaved to {aa_output_path}")
