import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("analysis/cid_comparison_destress_vs_aa_composition.csv")

output_path = "analysis/cid_destress_vs_aa_composition.png"

COLOR_BELOW = "#2a78d6"   # DE-STRESS better (diff < 0)
COLOR_ABOVE = "#eb6834"   # AA comp better (diff > 0)

diff = df["cid_diff_destress_minus_aa"].values   # negative = DE-STRESS closer
n_below = (diff < 0).sum()
n_above = (diff > 0).sum()

rng = np.random.default_rng(42)

fig, ax = plt.subplots(figsize=(5, 5))

# violin
parts = ax.violinplot(diff, positions=[0], widths=0.5,
                      showmeans=False, showmedians=False, showextrema=False)
for pc in parts["bodies"]:
    pc.set_facecolor("#2a78d6")
    pc.set_alpha(0.2)
    pc.set_edgecolor("#2a78d6")
    pc.set_linewidth(1)

# jittered points coloured by direction
jitter = rng.uniform(-0.08, 0.08, size=len(diff))
colors = [COLOR_BELOW if d < 0 else COLOR_ABOVE for d in diff]
ax.scatter(jitter, diff, c=colors, alpha=0.75, s=32, linewidths=0, zorder=3)

# mean + 1 SD
mean = diff.mean()
sd   = diff.std(ddof=1)
ax.errorbar(0, mean, yerr=sd,
            fmt="o", color="white",
            markersize=10, markeredgewidth=2, markeredgecolor="#2a78d6",
            elinewidth=2.5, capsize=8, capthick=2.5, ecolor="#2a78d6",
            zorder=5)

# zero line
ax.axhline(0, color="#333333", linewidth=1.2, linestyle="--", zorder=2)
ax.text(0.27, 0.004, "No difference", color="#555555", fontsize=8, va="bottom")

# annotations
ax.text(0.27,  diff.min() - 0.002, f"DE-STRESS closer  ({n_below}/{len(diff)})",
        color=COLOR_BELOW, fontsize=8, va="top")
ax.text(0.27,  diff.max() + 0.002, f"AA comp closer  ({n_above}/{len(diff)})",
        color=COLOR_ABOVE, fontsize=8, va="bottom")

pad = 0.015
ax.set_ylim(diff.min() - pad * 4, diff.max() + pad * 4)
ax.set_xlim(-0.4, 0.4)
ax.set_xticks([])
ax.set_ylabel("CID difference  (DE-STRESS minus AA composition)", fontsize=10)
ax.yaxis.grid(True, linewidth=0.6, color="#e0e0e0", zorder=0)
ax.set_axisbelow(True)
ax.spines[["top", "right", "bottom"]].set_visible(False)

ax.set_title(
    "DE-STRESS vs AA composition - paired CID difference\n"
    "Each point: one linkage × distance × scaler combination (n=39)",
    fontsize=10, pad=10,
)

plt.tight_layout()
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved to {output_path}")
