import pandas as pd

output_path = "analysis/data_exploration/destress_features_by_dataset.csv"

# Use standard scaling as representative - feature selection is consistent across scaling methods
# mode "cols": final features are the columns of the scaled CSV
# mode "rows": final features are listed row-wise (scFv feature selection outputs)
analyses = [
    ("PCA - AFDB full",                                    "AFDB pLDDT > 70",                 "data/processed_data/af2/standard/processed_destress_data_scaled.csv",                        "cols"),
    ("PCA - AFDB non-redundant (avg by organism)",                           "AFDB pLDDT > 70 & non-redundant", "data/processed_data/af2/standard/processed_destress_data_scaled_nonredundant.csv", "cols"),
    ("PCA - AFDB non-redundant (avg by organism and subcellular location)",  "AFDB pLDDT > 70 & non-redundant", "data/processed_data/af2/standard/processed_destress_data_scaled_nonredundant.csv", "cols"),
    ("PCA - PDB",                                          "PDB",                             "data/processed_data/pdb/standard/processed_destress_data_scaled.csv",                        "cols"),
    ("Hierarchical clustering - all organisms",            "AFDB pLDDT > 70 & non-redundant", "data/processed_data/af2/standard/processed_destress_data_scaled_nonredundant.csv",           "cols"),
    ("Hierarchical clustering - eukaryotes only",          "AFDB pLDDT > 70 & non-redundant", "data/processed_data/af2/standard/processed_destress_data_scaled_nonredundant.csv",           "cols"),
    ("Hierarchical clustering - 11 AFDB clusters",         "AFDB pLDDT > 70 (per cluster subset)", "data/processed_data/af2/standard/processed_destress_data_scaled.csv",                 "cols"),
    ("scFv classification - with composition (RF)",        "scFv Fleishman dataset",          "antibodyproduction/feature_selection/standard/comp/selected_features_rf.csv",                "rows"),
    ("scFv classification - with composition (MI)",        "scFv Fleishman dataset",          "antibodyproduction/feature_selection/standard/comp/selected_features_mi.csv",                "rows"),
    ("scFv classification - no composition (RF)",          "scFv Fleishman dataset",          "antibodyproduction/feature_selection/standard/no_comp/selected_features_rf.csv",             "rows"),
    ("scFv classification - no composition (MI)",          "scFv Fleishman dataset",          "antibodyproduction/feature_selection/standard/no_comp/selected_features_mi.csv",             "rows"),
]

rows = []
for analysis_name, dataset_name, path, mode in analyses:
    if mode == "cols":
        features = list(pd.read_csv(path, nrows=0).columns)
    else:
        features = list(pd.read_csv(path).iloc[:, 0])
    rows.append({
        "Analysis": analysis_name,
        "Dataset": dataset_name,
        "Features": ", ".join(sorted(features)),
        "N features": len(features),
    })

result = pd.DataFrame(rows)
result.to_csv(output_path, index=False)
print(result[["Analysis", "Dataset", "N features", "Features"]].to_string(index=False))
print(f"\nSaved to {output_path}")
