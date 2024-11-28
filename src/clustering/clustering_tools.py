# 0. Importing packages---------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import decomposition
from scipy.cluster.hierarchy import dendrogram


# 1. Defining helper functions--------------------------------------------


def adj_rand_ind_wssd_plot(data, title, file_name, output_path):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    sns.lineplot(
        data=data,
        x="n_clusters",
        y="weighted_ssd",
        errorbar="sd",
        legend=True,
        ax=ax1,
        color="blue",
        label="Weighted Sum of Square Distances",
    )
    sns.lineplot(
        data=data,
        x="n_clusters",
        y="adj_rand_score",
        errorbar="sd",
        legend=True,
        ax=ax2,
        color="red",
        label="Adjusted Rand Index",
    )
    ax1.set_xlabel("Number of clusters")
    ax1.set_ylabel("Weighted Sum of Square Distances")
    ax2.set_ylabel("Adjusted Rand Index")
    x_ticks = range(data["n_clusters"].min(), data["n_clusters"].max() + 1, 2)
    ax1.set_xticks(x_ticks)
    ax2.set_ylim([0, 1.2])
    plt.title(title)
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles = handles1 + handles2
    labels = labels1 + labels2
    plt.legend(handles, labels)
    # sns.move_legend(
    #     "upper left",
    #     bbox_to_anchor=(1, 1),
    #     frameon=False,
    # )
    plt.savefig(
        output_path + file_name + ".png",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close()


def adj_rand_ind_plot(data, title, file_name, output_path):
    # Setting theme for plots
    sns.set_style("whitegrid")
    sns.lineplot(
        data=data,
        x="n_clusters",
        y="adj_rand_score",
        # hue=hue,
        errorbar="sd",
        legend=True,
        lw=2,
    )
    plt.xlabel("Number of clusters", fontsize=16)
    plt.ylabel("Adjusted Rand Index", fontsize=16)
    # x_ticks = range(data["n_clusters"].min(), data["n_clusters"].max() + 1, 2)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlim(2, 20)
    plt.ylim([0, 1])
    plt.title(title)
    plt.savefig(
        output_path + file_name + ".png",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close()


def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # hierarchy.set_link_color_palette(sns.color_palette("colorblind", as_cmap=True))

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)


def get_newick(node, parent_dist, leaf_names, newick="") -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(
            node.get_right(), node.dist, leaf_names, newick=",%s" % (newick)
        )
        newick = "(%s" % (newick)
        return newick
