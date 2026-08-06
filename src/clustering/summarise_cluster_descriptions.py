import pandas as pd

labels_df = pd.read_csv("data/processed_data/af2/labels.csv")

org_counts = (
    labels_df.groupby("cluster_representative")["organism_scientific_name"]
    .nunique()
    .reset_index()
    .rename(columns={"organism_scientific_name": "n_orgs"})
)
af2db_cluster_list_df = org_counts[org_counts["n_orgs"] >= 40].reset_index(drop=True)

def get_top_desc(cluster_id):
    descs = labels_df["uniprot_description"][labels_df["cluster_representative"] == cluster_id]
    return descs.value_counts().index[0]

af2db_cluster_list_df["top_description"] = af2db_cluster_list_df["cluster_representative"].apply(get_top_desc)

exclude_pattern = "Uncharacterized protein|domain-containing protein"
af2db_cluster_list_df = af2db_cluster_list_df[
    ~af2db_cluster_list_df["top_description"].str.contains(exclude_pattern, na=False)
]

rows = []
for _, row in af2db_cluster_list_df.iterrows():
    cluster_mask = labels_df["cluster_representative"] == row["cluster_representative"]
    total = cluster_mask.sum()
    counts = labels_df.loc[cluster_mask, "uniprot_description"].value_counts()
    for desc, n in counts.items():
        is_variant = row["top_description"].lower() in desc.lower()
        rows.append({
            "cluster_representative": row["cluster_representative"],
            "cluster_name": row["top_description"],
            "uniprot_description": desc,
            "count": n,
            "total_in_cluster": total,
            "fraction": round(n / total, 3),
            "name_group": "top description variant" if is_variant else "other",
        })

out_df = pd.DataFrame(rows)
out_path = "analysis/hier_clustering_avg_by_org_single_af2db_cluster_filtered/cluster_description_counts.csv"
out_df.to_csv(out_path, index=False)
print(f"Saved to {out_path} ({len(out_df)} rows)")
