import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

PCA_CSV = "analysis/pca_single_af2db_cluster_filtered/af2/standard/A0A0G9LKG2/pca_transformed_data.csv"
OUTPUT_DIR = "analysis/pca_single_af2db_cluster_filtered/af2/standard/A0A0G9LKG2/"

df = pd.read_csv(PCA_CSV)

groups = {"Eukaryotes": df[df["organism_group2"] == "Eukaryotes"],
          "Prokaryotes": df[df["organism_group2"] == "Prokaryotes"]}

colors = {"Eukaryotes": "#0173b2", "Prokaryotes": "#de8f05"}

print("PC1 vs mean pLDDT within groups:\n")
for group_name, gdf in groups.items():
    pearson_r, pearson_p = stats.pearsonr(gdf["Mean_PLDDT"], gdf["dim0"])
    spearman_r, spearman_p = stats.spearmanr(gdf["Mean_PLDDT"], gdf["dim0"])
    print(f"  {group_name} (n={len(gdf)}, pLDDT range {gdf['Mean_PLDDT'].min():.1f}–{gdf['Mean_PLDDT'].max():.1f}):")
    print(f"    Pearson  r = {pearson_r:.3f}  p = {pearson_p:.3e}")
    print(f"    Spearman r = {spearman_r:.3f}  p = {spearman_p:.3e}")
    print()

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

for ax, (group_name, gdf) in zip(axes, groups.items()):
    pearson_r, pearson_p = stats.pearsonr(gdf["Mean_PLDDT"], gdf["dim0"])
    spearman_r, spearman_p = stats.spearmanr(gdf["Mean_PLDDT"], gdf["dim0"])

    ax.scatter(gdf["Mean_PLDDT"], gdf["dim0"],
               alpha=0.7, s=30, color=colors[group_name], label=group_name)

    m, b = np.polyfit(gdf["Mean_PLDDT"], gdf["dim0"], 1)
    x_line = np.linspace(gdf["Mean_PLDDT"].min(), gdf["Mean_PLDDT"].max(), 100)
    ax.plot(x_line, m * x_line + b, color="black", linewidth=1.2, linestyle="--")

    ax.set_xlabel("Mean pLDDT", fontsize=11)
    ax.set_ylabel("PC1", fontsize=11)
    ax.set_title(
        f"{group_name} (n={len(gdf)})\n"
        f"Pearson r = {pearson_r:.2f} (p = {pearson_p:.2e})\n"
        f"Spearman ρ = {spearman_r:.2f} (p = {spearman_p:.2e})",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.suptitle("SOD cluster: mean pLDDT vs PC1 within eukaryotes and prokaryotes", fontsize=11)
plt.tight_layout()

out_png = OUTPUT_DIR + "plddt_pca_correlation_by_kingdom.png"
plt.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"Saved to {out_png}")
