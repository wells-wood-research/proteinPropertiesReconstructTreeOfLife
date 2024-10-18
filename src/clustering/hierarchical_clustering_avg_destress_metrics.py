# 0. Importing packages---------------------------------------------------------

from scipy.cluster import hierarchy
from clustering_tools import *
import ete3
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import dendrogram
import dendropy

# 1. Defining variables---------------------------------------------------------

# Defining the scaling methods list
scaling_method_list = ["standard", "robust", "minmax"]

# Creating a color palette
palette = sns.color_palette(
    ["#0173b2", "#d55e00", "#029e73", "#cc78bc", "#808080", "#f0e442"], 6
)

# Defining data path
data_path = "data/processed_data/af2/"
raw_data_path = "data/raw_data/"

# Defining output data path
output_path = "analysis/hier_clustering_avg_by_org/af2/"

# Defining a dictionary of labels
label_dict = {
    "organism_group": "Organism Group 1",
    "organism_group2": "Organism Group 2",
}

# Defining the different linkage metrics for hieracrhcial clustering
linkage_list = ["single", "average", "complete", "ward"]

# Defining the different distance metrics
distance_metric_list = ["euclidean", "cityblock", "cosine", "correlation"]

# 2. Looping through the different scaling methods--------------------------------------------------

for scaling_method in scaling_method_list:

    # Defining paths for the data and scaling method used
    data_path_scaled = data_path + scaling_method + "/"

    # Defining the path for processed AF2 DE-STRESS data
    processed_destress_data_path = (
        data_path_scaled + "processed_destress_data_scaled.csv"
    )

    # Defining file paths for labels
    labels_df_path = data_path + "labels.csv"

    # 3. Reading in data------------------------------------------------------------------------

    # Defining the path for the processed AF2 DE-STRESS data
    processed_destress_data = pd.read_csv(processed_destress_data_path)

    # Reading in labels
    labels_df = pd.read_csv(labels_df_path)

    # Joining on the organism label data
    processed_destress_data_joined = pd.concat(
        [
            processed_destress_data,
            labels_df[
                ["organism_scientific_name", "organism_group", "organism_group2"]
            ],
        ],
        axis=1,
    )

    # # Filtering for euk
    # processed_destress_data_joined = processed_destress_data_joined[
    #     processed_destress_data_joined["organism_group2"] == "Eukaryotes"
    # ].reset_index(drop=True)

    # Average each principal component grouped by organism
    processed_destress_data_avg = processed_destress_data_joined.groupby(
        ["organism_scientific_name", "organism_group", "organism_group2"],
        as_index=False,
    )[processed_destress_data.columns.to_list()].mean()

    processed_destress_data_avg.update(
        processed_destress_data_avg[["organism_scientific_name"]].map("'{}'".format)
    )

    # Extracting labels
    organism_group_labels = processed_destress_data_avg["organism_group"].to_list()
    organism_group2_labels = processed_destress_data_avg["organism_group2"].to_list()
    organism_labels = processed_destress_data_avg["organism_scientific_name"].to_list()

    # Extracting labels
    labels = processed_destress_data_avg[
        ["organism_scientific_name", "organism_group", "organism_group2"]
    ]

    # Removing these labels from destress data
    processed_destress_data_avg.drop(
        ["organism_scientific_name", "organism_group", "organism_group2"],
        inplace=True,
        axis=1,
    )

    # 4. Running different hierarchical clustering--------------------------------------------

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
