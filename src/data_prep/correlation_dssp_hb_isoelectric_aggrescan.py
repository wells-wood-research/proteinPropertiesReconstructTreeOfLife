import pandas as pd
from scipy import stats

HB_COLS = ["rosetta_hbond_sr_bb", "rosetta_hbond_lr_bb"]
AGG_COLS = ["aggrescan3d_avg_value", "aggrescan3d_total_value", "aggrescan3d_min_value", "aggrescan3d_max_value"]
VDW_COLS = ["rosetta_fa_atr", "rosetta_fa_rep"]

labels = pd.read_csv("data/processed_data/af2/labels_nonredundant.csv").drop(columns=["aggrescan3d_avg_value", "aggrescan3d_avg_bin"], errors="ignore")
raw = pd.read_csv("data/raw_data/destress_data_af2.csv", usecols=["design_name", "num_residues"] + HB_COLS + AGG_COLS + VDW_COLS)

df = labels.merge(raw, on="design_name", how="inner")

# Normalise HB and VdW energies by sequence length to remove protein size effect
for col in HB_COLS + VDW_COLS:
    df[col + "_per_res"] = df[col] / df["num_residues"]

hb_per_res = [c + "_per_res" for c in HB_COLS]
vdw_per_res = [c + "_per_res" for c in VDW_COLS]

print("=" * 70)
print("DSSP bin vs Rosetta HB energies (per residue)")
print("Kruskal-Wallis test across DSSP bins")
print("=" * 70)

for col in hb_per_res:
    groups = [grp[col].dropna().values for _, grp in df.groupby("dssp_bin")]
    stat, p = stats.kruskal(*groups)
    print(f"\n{col}")
    print(f"  Kruskal-Wallis H={stat:.2f}, p={p:.2e}")
    medians = df.groupby("dssp_bin")[col].median().sort_values()
    for bin_label, median in medians.items():
        n = df[df["dssp_bin"] == bin_label][col].notna().sum()
        print(f"  {bin_label:<20} median={median:.4f}  (n={n})")

print()
print("=" * 70)
print("Isoelectric point vs Aggrescan3D metrics")
print("Spearman correlation (AF2 non-redundant)")
print("=" * 70)

agg_cols_available = [c for c in AGG_COLS if c in df.columns]
for col in agg_cols_available:
    subset = df[["isoelectric_point", col]].dropna()
    r, p = stats.spearmanr(subset["isoelectric_point"], subset[col])
    print(f"\nisoelectric_point vs {col}")
    print(f"  Spearman r={r:.3f}, p={p:.2e}  (n={len(subset)})")

print()
print("=" * 70)
print("Rosetta VdW energy (per residue) vs packing density")
print("Spearman correlation (AF2 non-redundant)")
print("=" * 70)

for col in vdw_per_res:
    subset = df[["packing_density", col]].dropna()
    r, p = stats.spearmanr(subset["packing_density"], subset[col])
    print(f"\npacking_density vs {col}")
    print(f"  Spearman r={r:.3f}, p={p:.2e}  (n={len(subset)})")
