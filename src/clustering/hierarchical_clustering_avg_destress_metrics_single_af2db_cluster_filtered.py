# 0. Importing packages---------------------------------------------------------

from scipy.cluster import hierarchy
from clustering_tools import *
import ete3
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import dendrogram
import dendropy
import os

# 1. Defining variables---------------------------------------------------------

# Defining the scaling methods list
scaling_method_list = ["standard", "robust", "minmax"]

# Creating a color palette
palette = sns.color_palette(
    ["#0173b2", "#d55e00", "#029e73", "#cc78bc", "#808080", "#f0e442"], 6
)

# Defining data path.
# Using the full processed dataset (untrimmed AF2 structures) across all clusters.
# Original single-cluster workflow (trimmed/aligned, A0A0G9LKG2 only):
# data_path = "data/processed_data_single_af2db_cluster/af2/"
data_path = "data/processed_data/af2/"
raw_data_path = "data/raw_data/"

# Defining output data path
output_path = "analysis/hier_clustering_avg_by_org_single_af2db_cluster_filtered/"

# Defining a path to the cluster list — FoldSeek clusters spanning >= 40 organisms
af2db_cluster_list_path = "data/processed_data/af2/uniq_org_counts_by_structural_cluster_gt_40.csv"

# Defining a dictionary of labels
label_dict = {
    "organism_group": "Organism Group 1",
    "organism_group2": "Organism Group 2",
}

# Defining the different linkage metrics for hieracrhcial clustering
linkage_list = ["single", "average", "complete", "ward"]

# Defining the different distance metrics
distance_metric_list = ["euclidean", "cityblock", "cosine", "correlation"]

# 2. Building the cluster list — done once, not per scaling method --------------------------

# Load labels (pLDDT-filtered processed data) to derive cluster list and descriptions.
labels_df = pd.read_csv(data_path + "labels.csv")

# Computing unique organism counts per cluster from the pLDDT-filtered processed data,
# so the >= 40 threshold is applied to the actual data being analysed (not pre-pLDDT counts).
org_counts = (
    labels_df.groupby("cluster_representative")["organism_scientific_name"]
    .nunique()
    .reset_index()
    .rename(columns={"organism_scientific_name": "n_orgs"})
)
af2db_cluster_list_df = org_counts[org_counts["n_orgs"] >= 40].reset_index(drop=True)

# Joining the most common UniProt description for each cluster from the labels data,
# so we can filter on protein identity rather than hardcoding cluster IDs.
labels_desc = labels_df[["cluster_representative", "uniprot_description"]]
top_desc = (
    labels_desc.groupby("cluster_representative")["uniprot_description"]
    .agg(lambda x: x.value_counts().index[0])
    .reset_index()
    .rename(columns={"uniprot_description": "top_description"})
)
af2db_cluster_list_df = af2db_cluster_list_df.merge(top_desc, on="cluster_representative", how="left")

# Excluding clusters that are not well-characterised single-function proteins.
# Two categories are excluded:
# 1. "Uncharacterized protein" - no functional interpretation is possible.
# 2. Domain-containing proteins (broad structural superfamilies, not true orthologous
#    families), making biological conclusions unreliable.
exclude_terms = ["Uncharacterized protein", "domain-containing protein"]
exclude_mask = af2db_cluster_list_df["top_description"].str.contains(
    "|".join(exclude_terms), na=False
)
af2db_cluster_list_df = af2db_cluster_list_df[~exclude_mask].reset_index(drop=True)
af2db_cluster_list = af2db_cluster_list_df["cluster_representative"].to_list()

print(f"Running hierarchical clustering for {len(af2db_cluster_list)} clusters:")
for c in af2db_cluster_list:
    desc = af2db_cluster_list_df.loc[af2db_cluster_list_df["cluster_representative"] == c, "top_description"].values[0]
    print(f"  {c}: {desc}")

# To run on all clusters including uncharacterised and domain-containing proteins:
# af2db_cluster_list = pd.read_csv(af2db_cluster_list_path)["cluster_representative"].to_list()

# 3. Looping through the different scaling methods--------------------------------------------------

for scaling_method in scaling_method_list:

    # Defining paths for the data and scaling method used
    data_path_scaled = data_path + scaling_method + "/"

    # Defining the path for processed AF2 DE-STRESS data
    processed_destress_data_path = (
        data_path_scaled + "processed_destress_data_scaled.csv"
    )

    # Defining file paths for labels
    labels_df_path = data_path + "labels.csv"

    # 4. Reading in data------------------------------------------------------------------------

    # Defining the path for the processed AF2 DE-STRESS data
    processed_destress_data = pd.read_csv(processed_destress_data_path)

    # Reading in labels
    labels_df = pd.read_csv(labels_df_path)

    # Joining on the organism label data (including cluster_representative for per-cluster filtering)
    processed_destress_data_joined = pd.concat(
        [
            processed_destress_data,
            labels_df[
                [
                    "organism_scientific_name",
                    "organism_group",
                    "organism_group2",
                    "subcellular_location",
                    "cluster_representative",
                ]
            ],
        ],
        axis=1,
    )

    # 5. Looping through each protein cluster ---------------------------------------------------

    for af2db_cluster in af2db_cluster_list:

        # Filter to structures belonging to this FoldSeek cluster.
        processed_destress_data_cluster = processed_destress_data_joined[
            processed_destress_data_joined["cluster_representative"] == af2db_cluster
        ].reset_index(drop=True)

        # Create per-cluster output directory: output_path/cluster/scaling_method/
        cluster_output_path = output_path + af2db_cluster + "/" + scaling_method + "/"
        os.makedirs(cluster_output_path, exist_ok=True)

        # Average each biophysical metric grouped by organism
        processed_destress_data_avg = processed_destress_data_cluster.groupby(
            ["organism_scientific_name"],
            as_index=False,
        )[processed_destress_data.columns.to_list()].mean()

        # Extracting organism name labels for dendrogram leaves
        organism_labels = processed_destress_data_avg["organism_scientific_name"].to_list()

        # Removing organism name column — clustering operates on numeric metrics only
        processed_destress_data_avg.drop(
            ["organism_scientific_name"],
            inplace=True,
            axis=1,
        )

        # 6. Running different hierarchical clustering--------------------------------------------

        for linkage in linkage_list:
            for distance_metric in distance_metric_list:

                if linkage == "ward" and distance_metric != "euclidean":
                    continue
                else:
                    linkage_matrix = hierarchy.linkage(
                        processed_destress_data_avg, method=linkage, metric=distance_metric
                    )

                    plt.figure(figsize=(9, 8))
                    dendrogram(
                        linkage_matrix,
                        truncate_mode=None,
                        labels=organism_labels,
                        orientation="left",
                        leaf_font_size=8,
                    )
                    plt.xticks(fontsize=10)
                    plt.savefig(
                        cluster_output_path
                        + "dendrogram_"
                        + scaling_method
                        + "_"
                        + "linkage-"
                        + linkage
                        + "_"
                        + "distance-"
                        + distance_metric
                        + ".png",
                        bbox_inches="tight",
                        dpi=600,
                    )
                    plt.close()

                    # Convert linkage matrix to dendrogram
                    dendro = dendrogram(linkage_matrix, no_plot=True)

                    # Convert the linkage matrix to a tree object
                    tree = hierarchy.to_tree(linkage_matrix, False)

                    # Convert the tree object to the Newick format
                    newick = get_newick(tree, tree.dist, organism_labels)

                    # Using ETE3
                    ete3_tree = ete3.Tree(newick)

                    # Save the tree in Newick format
                    ete3_tree.write(
                        format=1,
                        outfile=cluster_output_path
                        + "tree_"
                        + scaling_method
                        + "_"
                        + "linkage-"
                        + linkage
                        + "_"
                        + "distance-"
                        + distance_metric
                        + ".nwk",
                    )
