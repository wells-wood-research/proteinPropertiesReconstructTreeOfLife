import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

PCA_CSV = "analysis/pca_single_af2db_cluster_filtered/af2/standard/A0A0G9LKG2/pca_transformed_data.csv"
OUTPUT_DIR = "analysis/pca_single_af2db_cluster_filtered/af2/standard/A0A0G9LKG2/"

df = pd.read_csv(PCA_CSV)

print(f"Structures: {len(df)}")
print(f"Mean pLDDT range: {df['Mean_PLDDT'].min():.1f} – {df['Mean_PLDDT'].max():.1f}")
print()

for pc_col, pc_label in [("dim0", "PC1"), ("dim1", "PC2")]:
    pearson_r, pearson_p = stats.pearsonr(df["Mean_PLDDT"], df[pc_col])
    spearman_r, spearman_p = stats.spearmanr(df["Mean_PLDDT"], df[pc_col])
    print(f"{pc_label} vs mean pLDDT:")
    print(f"  Pearson  r = {pearson_r:.3f}  p = {pearson_p:.3e}")
    print(f"  Spearman r = {spearman_r:.3f}  p = {spearman_p:.3e}")
    print()

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

for ax, pc_col, pc_label in zip(axes, ["dim0", "dim1"], ["PC1", "PC2"]):
    pearson_r, pearson_p = stats.pearsonr(df["Mean_PLDDT"], df[pc_col])
    spearman_r, spearman_p = stats.spearmanr(df["Mean_PLDDT"], df[pc_col])

    ax.scatter(df["Mean_PLDDT"], df[pc_col], alpha=0.6, s=20, color="#0173b2")

    m, b = np.polyfit(df["Mean_PLDDT"], df[pc_col], 1)
    x_line = np.linspace(df["Mean_PLDDT"].min(), df["Mean_PLDDT"].max(), 100)
    ax.plot(x_line, m * x_line + b, color="#de8f05", linewidth=1.5)

    ax.set_xlabel("Mean pLDDT", fontsize=11)
    ax.set_ylabel(pc_label, fontsize=11)
    ax.set_title(
        f"r = {pearson_r:.2f} (p = {pearson_p:.2e})\n"
        f"Spearman ρ = {spearman_r:.2f} (p = {spearman_p:.2e})",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.suptitle("SOD cluster: mean pLDDT vs PC scores (standard scaling)", fontsize=11)
plt.tight_layout()

out_png = OUTPUT_DIR + "plddt_pca_correlation.png"
plt.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"Saved plot to {out_png}")
