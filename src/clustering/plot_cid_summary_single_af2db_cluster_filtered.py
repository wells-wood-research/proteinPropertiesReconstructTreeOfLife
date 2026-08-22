import textwrap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set_style("whitegrid")

summary_csv = "analysis/hier_clustering_avg_by_org_single_af2db_cluster_filtered/tree_distances_all_clusters.csv"
output_path = "analysis/hier_clustering_avg_by_org_single_af2db_cluster_filtered/"

df = pd.read_csv(summary_csv)

# Per-cluster summary
stats = (
    df.groupby("description")
    .agg(
        mean_cid=("tree_dist", "mean"),
        sd_cid_obs=("tree_dist", "std"),
        expected_cid=("expected_cid", "mean"),
        expected_sd=("sd_cid", "mean"),
    )
    .reset_index()
    .sort_values("mean_cid", ascending=True)
)

stats_path = output_path + "cid_summary_stats.csv"
stats.to_csv(stats_path, index=False)
print(f"Saved cluster stats to {stats_path}")

fig, ax = plt.subplots(figsize=(11, 5))

y_positions = np.arange(len(stats))
descriptions = stats["description"].tolist()

# Individual CID values as jittered background points
rng = np.random.default_rng(42)
for i, desc in enumerate(descriptions):
    vals = df.loc[df["description"] == desc, "tree_dist"].values
    jitter = rng.uniform(-0.2, 0.2, size=len(vals))
    ax.scatter(vals, i + jitter, color="#0173b2", alpha=0.25, s=10, zorder=3)

# Per-cluster expected CID: diamond marker with ±1 SD caps
for i, (_, row) in enumerate(stats.iterrows()):
    ax.errorbar(
        row["expected_cid"], i,
        xerr=row["expected_sd"],
        fmt="D",
        color="#999999",
        markersize=5,
        linewidth=1.0,
        capsize=3,
        zorder=2,
    )

# Observed mean ± SD error bars and mean points
ax.errorbar(
    stats["mean_cid"],
    y_positions,
    xerr=stats["sd_cid_obs"],
    fmt="o",
    color="#0173b2",
    markersize=6,
    linewidth=1.2,
    capsize=3,
    zorder=4,
)

wrapped = ["\n".join(textwrap.wrap(d, width=35)) for d in descriptions]
ax.set_yticks(y_positions)
ax.set_yticklabels(wrapped, fontsize=9)
ax.set_xlabel("Clustering Information Distance (CID)", fontsize=11)
ax.set_xlim(0.45, 0.95)
ax.set_xticks(np.arange(0.5, 0.95, 0.1))

# Legend
jitter_handle = plt.Line2D([0], [0], marker="o", linestyle="none",
                            markerfacecolor="#0173b2", markeredgecolor="none",
                            markersize=5, alpha=0.4,
                            label="Observed CID values (DE-STRESS tree vs. NCBI reference)")
obs_handle = plt.Line2D([0], [0], marker="o", color="#0173b2", linestyle="-",
                         markersize=6, label="Mean ± SD of observed CID values")
exp_handle = plt.Line2D([0], [0], marker="D", color="#999999", linestyle="-",
                         markersize=5, label="Mean ± SD of random tree CID values")
ax.legend(handles=[jitter_handle, obs_handle, exp_handle], fontsize=8,
          loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, borderaxespad=0)


plt.tight_layout()

out_file = output_path + "cid_summary_plot.png"
plt.savefig(out_file, dpi=300, bbox_inches="tight")
print(f"Saved CID summary plot to {out_file}")

out_pdf = output_path + "cid_summary_plot.pdf"
plt.savefig(out_pdf, bbox_inches="tight")
print(f"Saved CID summary plot to {out_pdf}")
