import pandas as pd
from scipy import stats

PCA_CSV   = "analysis/pca_single_af2db_cluster_filtered/af2/standard/A0A0G9LKG2/pca_transformed_data.csv"
FEAT_CSV  = "data/processed_data/af2/standard/processed_destress_data_scaled.csv"
LABEL_CSV = "data/processed_data/af2/labels.csv"

pca_df   = pd.read_csv(PCA_CSV)[["design_name", "dim0", "Mean_PLDDT", "organism_group2"]]
feat_df  = pd.read_csv(FEAT_CSV)
label_df = pd.read_csv(LABEL_CSV)[["design_name", "cluster_representative"]]

# Row-align features and labels (both are full dataset, same order)
combined = pd.concat([label_df, feat_df], axis=1)

# Filter to SOD cluster then merge on design_name to get PC scores
sod = combined[combined["cluster_representative"] == "A0A0G9LKG2"].copy()
merged = sod.merge(pca_df, on="design_name")

euk = merged[merged["organism_group2"] == "Eukaryotes"]

feature_cols = feat_df.columns.tolist()

print(f"Eukaryotes n={len(euk)}\n")
print(f"{'Feature':<40} {'Pearson r':>10} {'p-value':>12}")
print("-" * 65)

results = []
for col in feature_cols:
    r, p = stats.pearsonr(euk[col], euk["dim0"])
    results.append((col, r, p))

results.sort(key=lambda x: abs(x[1]), reverse=True)
for col, r, p in results:
    print(f"{col:<40} {r:>10.3f} {p:>12.3e}")

print()
plddt_r, plddt_p = stats.pearsonr(euk["Mean_PLDDT"], euk["dim0"])
print(f"{'Mean_PLDDT (not a PCA feature)':<40} {plddt_r:>10.3f} {plddt_p:>12.3e}")
n_stronger = sum(1 for _, r, _ in results if abs(r) > abs(plddt_r))
print(f"\nFeatures with stronger PC1 correlation than pLDDT: {n_stronger}/{len(results)}")
