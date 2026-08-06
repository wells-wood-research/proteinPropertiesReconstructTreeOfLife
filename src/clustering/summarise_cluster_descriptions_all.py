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

def get_inclusion(top_desc):
    if "Uncharacterized protein" in top_desc:
        return "excluded", "Top annotation is uncharacterized protein — no functional identity"
    if "domain-containing protein" in top_desc:
        return "excluded", "Domain-level annotation only — structural fold name, not a specific protein"
    return "included", "Specific named protein"

af2db_cluster_list_df[["included", "inclusion_reason"]] = pd.DataFrame(
    af2db_cluster_list_df["top_description"].apply(get_inclusion).tolist(),
    index=af2db_cluster_list_df.index,
)

rows = []
for _, row in af2db_cluster_list_df.iterrows():
    cluster_mask = labels_df["cluster_representative"] == row["cluster_representative"]
    total = cluster_mask.sum()
    counts = labels_df.loc[cluster_mask, "uniprot_description"].value_counts()
    n_variant = sum(
        n for desc, n in counts.items()
        if row["top_description"].lower() in desc.lower()
    )
    value_counts_str = "; ".join(f"{desc}: {n}" for desc, n in counts.items())
    cluster_data = labels_df.loc[cluster_mask]
    organisms_str = "; ".join(sorted(cluster_data["organism_scientific_name"].dropna().unique()))
    kingdoms_str = "; ".join(sorted(cluster_data["organism_group"].dropna().unique()))
    rows.append({
        "cluster_representative": row["cluster_representative"],
        "cluster_name": row["top_description"],
        "n_orgs": row["n_orgs"],
        "total_proteins": total,
        "pct_matching_top_description": round(100 * n_variant / total, 1),
        "included": row["included"],
        "inclusion_reason": row["inclusion_reason"],
        "description_value_counts": value_counts_str,
        "organisms": organisms_str,
        "kingdoms": kingdoms_str,
    })

out_df = pd.DataFrame(rows).sort_values("pct_matching_top_description", ascending=False)
out_path = "analysis/hier_clustering_avg_by_org_single_af2db_cluster_filtered/cluster_summary_all22.csv"
out_df.to_csv(out_path, index=False)
print(f"Saved to {out_path}")
print(out_df.to_string(index=False))
