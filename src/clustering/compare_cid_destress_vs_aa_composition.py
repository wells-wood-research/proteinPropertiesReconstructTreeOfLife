import pandas as pd

destress = pd.read_csv("analysis/hier_clustering_avg_by_org/af2/tree_distances.csv")
aa_comp  = pd.read_csv("analysis/hier_clustering_avg_by_org_aa_composition/af2/tree_distances.csv")

output_path = "analysis/cid_comparison_destress_vs_aa_composition.csv"

# R converts hyphens to dots in column names, so both CSVs use dots
# Format: tree_<scaler>_linkage.<method>_distance.<metric>.nwk
def parse_name(name):
    name = name.replace(".nwk", "")
    parts = name.split("_")
    scaler   = parts[1]
    linkage  = parts[2].split(".")[1]
    distance = parts[3].split(".")[1]
    return scaler, linkage, distance

d = destress.T.reset_index()
d.columns = ["tree", "cid_destress"]
d[["scaler", "linkage", "distance"]] = d["tree"].apply(
    lambda x: pd.Series(parse_name(x))
)

a = aa_comp.T.reset_index()
a.columns = ["tree", "cid_aa_comp"]
a[["scaler", "linkage", "distance"]] = a["tree"].apply(
    lambda x: pd.Series(parse_name(x))
)

m = d.merge(a, on=["scaler", "linkage", "distance"])
m["cid_destress"] = m["cid_destress"].round(4)
m["cid_aa_comp"]  = m["cid_aa_comp"].round(4)
m["cid_diff_destress_minus_aa"] = (m["cid_destress"] - m["cid_aa_comp"]).round(4)
m["better_feature_set"] = m["cid_diff_destress_minus_aa"].apply(
    lambda x: "DE-STRESS" if x < 0 else "AA composition"
)

out = m[
    ["scaler", "linkage", "distance",
     "cid_destress", "cid_aa_comp",
     "cid_diff_destress_minus_aa", "better_feature_set"]
].sort_values(["linkage", "scaler", "distance"]).reset_index(drop=True)

out.to_csv(output_path, index=False)
print(f"Saved to {output_path}")
print(f"Rows: {len(out)}")
print()
print(out.to_string(index=False))
print()

random_mean = 0.8862
random_sd   = 0.0156
print(f"Random baseline (48 tips): mean={random_mean}, sd={random_sd}")
print()

summary = pd.DataFrame({
    "feature_set": ["DE-STRESS", "AA composition", "Random baseline"],
    "min_cid":  [out["cid_destress"].min(),          out["cid_aa_comp"].min(),  random_mean],
    "mean_cid": [round(out["cid_destress"].mean(), 4), round(out["cid_aa_comp"].mean(), 4), random_mean],
    "max_cid":  [out["cid_destress"].max(),           out["cid_aa_comp"].max(),  random_mean],
})
print(summary.to_string(index=False))
print()

n_total = len(out)
n_destress_better = (out["cid_diff_destress_minus_aa"] < 0).sum()
print(f"DE-STRESS outperforms AA composition: {n_destress_better}/{n_total} combinations")
