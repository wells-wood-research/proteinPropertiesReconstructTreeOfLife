# Hierarchical clustering baseline using amino acid composition (20 standard AAs).
# Reads scaled AA composition data produced by data_prep_aa_composition.py,
# averages per organism, and runs all linkage x distance combinations.
# Outputs dendrograms and Newick .nwk files to
# analysis/hier_clustering_avg_by_org_aa_composition/af2/.

# 0. Importing packages---------------------------------------------------------

from scipy.cluster import hierarchy
from clustering_tools import *
import ete3
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import dendrogram
import dendropy

# 1. Defining variables---------------------------------------------------------

scaling_method_list = ["standard", "robust", "minmax"]

data_path = "data/processed_data_aa_composition/af2/"

output_path = "analysis/hier_clustering_avg_by_org_aa_composition/af2/"

linkage_list = ["single", "average", "complete", "ward"]

distance_metric_list = ["euclidean", "cityblock", "cosine", "correlation"]

# 2. Looping through the different scaling methods------------------------------

for scaling_method in scaling_method_list:

    data_path_scaled = data_path + scaling_method + "/"

    scaled_data_path = data_path_scaled + "aa_composition_scaled_nonredundant.csv"
    labels_df_path = data_path + "labels_nonredundant.csv"

    # 3. Reading in data--------------------------------------------------------

    scaled_data = pd.read_csv(scaled_data_path)
    labels_df = pd.read_csv(labels_df_path)

    data_joined = pd.concat(
        [
            scaled_data,
            labels_df[["organism_scientific_name", "organism_group", "organism_group2"]],
        ],
        axis=1,
    )

    # Average per organism
    feature_cols = scaled_data.columns.to_list()

    data_avg = data_joined.groupby(
        ["organism_scientific_name"],
        as_index=False,
    )[feature_cols].mean()

    organism_labels = data_avg["organism_scientific_name"].to_list()

    data_avg.drop(["organism_scientific_name"], inplace=True, axis=1)

    # 4. Running hierarchical clustering----------------------------------------

    for linkage in linkage_list:
        for distance_metric in distance_metric_list:

            if linkage == "ward" and distance_metric != "euclidean":
                continue

            linkage_matrix = hierarchy.linkage(
                data_avg, method=linkage, metric=distance_metric
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
                output_path
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

            # Convert linkage matrix to Newick tree
            dendrogram(linkage_matrix, no_plot=True)

            tree = hierarchy.to_tree(linkage_matrix, False)

            # Single-quote labels so spaces are preserved when ete3 reads the Newick
            quoted_labels = ["'" + label + "'" for label in organism_labels]
            newick = get_newick(tree, tree.dist, quoted_labels)

            ete3_tree = ete3.Tree(newick)

            ete3_tree.write(
                format=1,
                outfile=output_path
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
