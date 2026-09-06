# Reads a tree_distances.csv produced by tree_dist_single_af2db_cluster_filtered.R,
# parses scaling method, linkage method, and distance metric from the file path,
# rounds CID to 2 d.p., sorts ascending by CID, and writes a formatted CSV.

import os
import re
import pandas as pd

input_path = (
    "analysis/hier_clustering_avg_by_org_single_af2db_cluster_filtered"
    "/A0A0G9LKG2/tree_distances.csv"
)

output_path = (
    "analysis/hier_clustering_avg_by_org_single_af2db_cluster_filtered"
    "/A0A0G9LKG2/tree_distances_formatted.csv"
)

df = pd.read_csv(input_path)

# Parse hyperparameters from file path
pattern = re.compile(
    r"/(standard|robust|minmax)/tree_(?:standard|robust|minmax)"
    r"_linkage-(\w+)_distance-(\w+)\.nwk"
)

rows = []
for _, row in df.iterrows():
    m = pattern.search(row["file"])
    if not m:
        raise ValueError(f"Could not parse: {row['file']}")
    scaling, linkage, distance = m.group(1), m.group(2), m.group(3)
    rows.append({
        "Scaling Method": scaling,
        "Linkage Method": linkage,
        "Distance Metric": distance,
        "Clustering Information Distance": round(row["tree_dist"], 2),
    })

out = (
    pd.DataFrame(rows)
    .sort_values("Clustering Information Distance")
    .reset_index(drop=True)
)

out.to_csv(output_path, index=False)
print(f"Written {len(out)} rows to {output_path}")
