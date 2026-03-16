# This script provides helper functions for the dimensionality reduction
# analysis.

# 0. Importing packages--------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
from sklearn import decomposition
from sklearn.ensemble import IsolationForest
from sklearn.metrics import pairwise_distances

# 1. Defining helper functions for Principal Component Analysis (PCA)---------------------------------


# Defining a function to calculate the variance explained by a
# number of different principal components
def pca_var_explained(data, n_components, file_name, output_path):
    # Performing PCA with the specified components
    pca = decomposition.PCA(n_components=n_components)
    pca.fit(data)

    # Calculating the variance explained
    var_explained = pca.explained_variance_ratio_

    # Calculating the cumulative sum
    var_explained_sum = np.cumsum(var_explained)

    # Calculating list of components
    components_list = range(1, n_components + 1, 1)

    # Creating dict
    var_explained_dict = {
        "n_components": components_list,
        "var_explained": var_explained,
        "var_explained_sum": var_explained_sum,
    }

    # Appending to the data frame
    var_explained_df = pd.DataFrame(var_explained_dict)

    # Saving as a csv file
    var_explained_df.to_csv(
        output_path + file_name + ".csv",
        index=False,
    )

    sns.set_style("whitegrid")

    # Plotting the data and saving
    sns.lineplot(
        x="n_components",
        y="var_explained_sum",
        data=var_explained_df,
    )
    # plt.title("""Cumulative variance explained by number of principal components""")
    plt.xlabel("Number of components", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel("Cumulative variance explained", fontsize=14)
    plt.savefig(output_path + file_name + ".png")
    plt.close()

    return var_explained_df


# Defining a script to perform Principal Component Analysis (PCA)
# for a data set and a specified number of principal components
def perform_pca(
    data, labels_df, n_components, output_path, file_path, components_file_path
):
    # Performing PCA
    pca_model = decomposition.PCA(n_components=n_components)
    pca_model.fit(data)

    # Saving contributions of the features to the principal components
    pca_feat_contr_to_cmpts = pd.DataFrame(
        np.round(abs(pca_model.components_), 4), columns=data.columns
    )

    pca_feat_contr_to_cmpts.to_csv(
        output_path + components_file_path + "feat_contr_to_cmpts.csv", index=True
    )

    # Defining the columns dict
    columns_dict = {}

    # Selecting the 10 largest contributers to each principal component
    for i in range(0, n_components):
        pca_components_contr = pca_feat_contr_to_cmpts.iloc[i].nlargest(
            10, keep="first"
        )
        pca_components_contr.to_csv(
            output_path + components_file_path + str(i) + "_contr.csv", index=True
        )

        columns_dict[i] = "dim" + str(i)

    # Transforming the data
    pca_transformed_data = pca_model.transform(data)

    # Converting to data frame and renaming columns
    pca_transformed_data = pd.DataFrame(pca_transformed_data).rename(
        columns=columns_dict
    )

    if labels_df is None:
        pca_transformed_data = pca_transformed_data

    else:
        # Adding the labels back
        pca_transformed_data = pd.concat([labels_df, pca_transformed_data], axis=1)

    # Outputting the transformed data
    pca_transformed_data.to_csv(
        output_path + file_path + ".csv",
        index=False,
    )

    return pca_transformed_data, pca_model


def plot_latent_space_2d(
    data,
    var_explained_data,
    x,
    y,
    axes_prefix,
    legend_title,
    hue,
    hue_order,
    alpha,
    s,
    palette,
    output_path,
    file_name,
):
    sns.color_palette("tab10")

    x_id = str(int(x[-1]) + 1)
    y_id = str(int(y[-1]) + 1)

    x_var_explained = var_explained_data["var_explained"][
        var_explained_data["n_components"] == 1
    ]
    y_var_explained = var_explained_data["var_explained"][
        var_explained_data["n_components"] == 2
    ]

    x_var_explained_formatted = np.round(x_var_explained.iloc[0], 2) * 100
    y_var_explained_formatted = np.round(y_var_explained.iloc[0], 2) * 100

    # PCA 2d scatter plot
    plot = sns.scatterplot(
        x=x,
        y=y,
        data=data,
        hue=hue,
        hue_order=hue_order,
        # style=style,
        alpha=alpha,
        palette=palette,
        s=s,
        legend=True,
        linewidth=0.2,
        edgecolor="black",
    )
    plt.xlabel(
        axes_prefix + x_id + " (" + str(np.int64(x_var_explained_formatted)) + "%)",
        fontsize=18,
    )
    plt.ylabel(
        axes_prefix + y_id + " (" + str(np.int64(y_var_explained_formatted)) + "%)",
        fontsize=18,
    )
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    # plt.xlim([-6, 11])
    # plt.ylim([-4, 8])
    # sns.move_legend(
    #     plot,
    #     "upper left",
    #     bbox_to_anchor=(1, 1),
    #     frameon=False,
    # )
    # sns.move_legend(
    #     plot,
    #     "lower center",
    #     bbox_to_anchor=(0.45, -0.35),
    #     frameon=False,
    #     ncols=6,
    #     fontsize=15,
    #     title_fontsize=15,
    # )
    sns.move_legend(
        plot,
        "lower center",
        bbox_to_anchor=(0.45, -0.40),
        frameon=True,
        ncols=4,
        title=legend_title,
        title_fontsize=14,
        fontsize=12,
    )
    plt.savefig(
        output_path + file_name + x_id + y_id + ".png",
        bbox_inches="tight",
        dpi=600,
    )
    plt.close()


def plot_pca_plotly(
    pca_data,
    x,
    y,
    z,
    color,
    hover_data,
    legend_title,
    opacity,
    size,
    output_path,
    file_name,
):
    fig = px.scatter_3d(
        pca_data,
        x=x,
        y=y,
        z=z,
        color=color,
        color_discrete_sequence=px.colors.qualitative.G10,
        hover_data=hover_data,
        opacity=opacity,
        labels={
            "dim0": "Principal Component 1",
            "dim1": "Principal Component 2",
            "dim2": "Principal Component 3",
        },
    )
    fig.update_traces(
        marker=dict(size=size, line=dict(width=0.8)),
        selector=dict(mode="markers"),
    )
    fig.update_layout(legend_title_text=legend_title)
    fig.write_html(output_path + file_name)


def distance_to_reference(data, dim_columns, feature, distance_metric, output_path):
    # Filtering
    pca_data_filt = data[dim_columns]

    # Computing distances
    distances = pairwise_distances(X=pca_data_filt, metric=distance_metric, n_jobs=1)

    # Converting to a data frame
    distances_df = pd.DataFrame(distances, columns=data[feature].to_list())

    distances_df = distances_df.round(decimals=4)

    # Outputting as a csv
    distances_df.to_csv(
        output_path + "pca" + str(len(dim_columns)) + "d_distances_all.csv", index=False
    )

    return distances_df


def spectral_plot(
    pca_data,
    group_var,
    hue_order,
    value_var_list,
    title,
    legend_title,
    output_path,
    file_name,
):
    # Changing format of data from wide to long
    pca_data_long = pca_data.melt(
        id_vars=[group_var],
        value_vars=value_var_list,
        var_name="dim_id",
        value_name="dim_value",
    )

    # Extracting the id of the pca dimension
    pca_data_long["dim_id"] = pca_data_long["dim_id"].str.replace("dim", "")

    sns.set_style("ticks")
    cmap = sns.color_palette("colorblind", as_cmap=True)

    # Plot the pca spectre for the different sequences
    plot = sns.lineplot(
        x="dim_id",
        y="dim_value",
        hue=group_var,
        hue_order=hue_order,
        # errorbar="sd",
        legend="full",
        data=pca_data_long,
        linewidth=4,
        palette=cmap,
        # alpha=0.9,
    )
    # ax.legend(
    #     loc="lower center",
    #     bbox_to_anchor=(0.5, -0.3),
    #     ncol=1,
    # )

    # sns.move_legend(
    #     plot,
    #     "upper left",
    #     bbox_to_anchor=(1, 1),
    #     frameon=False,
    #     fontsize=14,
    #     title=legend_title,
    # )
    # plot.get_legend().set_title(legend_title)
    # plt.legend(title=legend_title, fontsize="20", title_fontsize="14")
    sns.move_legend(
        plot,
        "lower center",
        bbox_to_anchor=(0.5, -0.38),
        frameon=False,
        ncols=2,
        fontsize=14,
        title=legend_title,
    )
    plt.xlabel("Principal Component ID", fontsize=14)
    plt.ylabel("Principal Component Value", fontsize=14)
    plt.title(title, fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.savefig(
        output_path + file_name + ".png",
        bbox_inches="tight",
        dpi=600,
    )
    plt.close()
