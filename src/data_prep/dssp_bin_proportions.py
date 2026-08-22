import pandas as pd

output_path = "analysis/data_exploration/dssp_bin_proportions.csv"

datasets = [
    ("AF2DB full",           "data/processed_data/af2/labels.csv"),
    ("AF2DB non-redundant",  "data/processed_data/af2/labels_nonredundant.csv"),
    ("PDB",                  "data/processed_data/pdb/labels.csv"),
]

# Row order for the output table
row_order = ["Mixed", "Alpha Helix", "Loop", "Beta Strand", "Bend", "Hbond Turn", "3 10 Helix"]

# Build a dict of {dssp_bin: {dataset: "X.X% (count)"}}
data = {}

for dataset_name, path in datasets:
    df = pd.read_csv(path)
    total = len(df)
    counts = df["dssp_bin"].value_counts()
    props = df["dssp_bin"].value_counts(normalize=True).mul(100)
    for bin_label in counts.index:
        if bin_label not in data:
            data[bin_label] = {}
        formatted = f"{props[bin_label]:.1f}% ({counts[bin_label]:,})"
        data[bin_label][dataset_name] = formatted

# Build rows in specified order, adding any unexpected categories at the end
ordered_bins = [b for b in row_order if b in data]
ordered_bins += [b for b in data if b not in row_order]

rows = []
for bin_label in ordered_bins:
    row = {"Secondary structure": bin_label}
    for dataset_name, _ in datasets:
        row[dataset_name] = data[bin_label].get(dataset_name, "0.0% (0)")
    rows.append(row)

result = pd.DataFrame(rows)
result.to_csv(output_path, index=False)
print(result.to_string(index=False))
print(f"\nSaved to {output_path}")
