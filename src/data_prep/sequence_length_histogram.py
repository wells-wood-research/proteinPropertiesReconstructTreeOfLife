import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

output_path = "analysis/data_exploration/sequence_length_histogram.png"

datasets = [
    ("AFDB pLDDT > 70",                 "data/processed_data/af2/labels.csv"),
    ("AFDB pLDDT > 70 & non-redundant", "data/processed_data/af2/labels_nonredundant.csv"),
    ("PDB",                             "data/processed_data/pdb/labels.csv"),
]

COLORS = ["#0173b2", "#d55e00", "#029e73"]

sns.set_style("whitegrid")

frames = []
for label, path in datasets:
    lengths = pd.read_csv(path, usecols=["full_sequence"])["full_sequence"].dropna().str.len()
    frames.append(pd.DataFrame({"length": lengths, "Dataset": label}))
df = pd.concat(frames, ignore_index=True)

palette = dict(zip([label for label, _ in datasets], COLORS))

fig, ax = plt.subplots(figsize=(9, 5))

sns.kdeplot(data=df, x="length", hue="Dataset", palette=palette,
            common_norm=False, common_grid=True, cut=0,
            fill=True, alpha=0.2, linewidth=2, ax=ax)

ax.set_xlabel("Sequence length (residues)", fontsize=11)
ax.set_ylabel("Density", fontsize=11)
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax.tick_params(labelsize=10)


plt.tight_layout()
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved to {output_path}")

for label, group in df.groupby("Dataset", sort=False):
    lengths = group["length"]
    print(f"{label}: n={len(lengths):,}  median={lengths.median():.0f}  mean={lengths.mean():.0f}  min={lengths.min()}  max={lengths.max():,}")
