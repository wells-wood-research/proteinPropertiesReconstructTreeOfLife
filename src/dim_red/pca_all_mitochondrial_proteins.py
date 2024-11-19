# This script runs Principal Component Analysis (PCA)
# across the different data sets and scaling methods.

# 0. Importing packages------------------------------------------------------------
from dim_red_tools import *
import os

# 1. Defining variables------------------------------------------------------------

# Defining the data set list
dataset_list = ["af2"]

# Defining the scaling methods list
# scaling_method_list = ["standard", "robust", "minmax"]
scaling_method_list = ["standard"]

# Defining number of principal components
n_components = 7

# Defining list of dim ids
dim_ids_list = []
for i in range(0, n_components):
    dim_ids_list.append("dim" + str(i))


# Defining hiver data for plotly
hover_data = ["organism_scientific_name", "dim0", "dim1", "organism_scientific_name"]

# Creating a color palette
palette = sns.color_palette(
    ["#0173b2", "#d55e00", "#029e73", "#cc78bc", "#808080", "#f0e442"], 6
)

# palette = sns.color_palette("tab10")

# Defining data path
data_path = "data/processed_data/"
raw_data_path = "data/raw_data/"

# Defining output data path
output_path = "analysis/pca_all_mitochondrial_proteins/"

# Defining a dictionary of labels
label_dict = {
    "organism_group_mito": "Organism Group 1",
    # "organism_group2_mito": "Organism Group 2",
}


# 2. Looping through the different data sets--------------------------------------------------

for dataset in dataset_list:

    # Defining the paths for the specified data set
    data_path_dataset = data_path + dataset + "/"
    output_path_dataset = output_path + dataset + "/"

    for scaling_method in scaling_method_list:

        # Defining paths for the data and scaling method used
        data_path_scaled = data_path_dataset + scaling_method + "/"
        output_path_scaled = output_path_dataset + scaling_method + "/"

        # Defining the path for processed AF2 DE-STRESS data
        processed_destress_data_path = (
            data_path_scaled + "processed_destress_data_scaled.csv"
        )

        # Defining file paths for labels
        labels_df_path = data_path_dataset + "labels.csv"

        # 3. Reading in data------------------------------------------------------------------------

        # Defining the path for the processed AF2 DE-STRESS data
        processed_destress_data = pd.read_csv(processed_destress_data_path)

        # Reading in labels
        labels_df = pd.read_csv(labels_df_path)

        # Filtering destress data
        processed_destress_data = processed_destress_data[
            (labels_df["Mean_PLDDT"] >= 70)
        ].reset_index(drop=True)

        # Filtering labels data
        labels_df = labels_df[(labels_df["Mean_PLDDT"] >= 70)].reset_index(drop=True)

        # Filtering for gene_encoding_type 'Mitochondrion'
        mitochondrion_df = processed_destress_data[
            labels_df["gene_encoding_type"] == "Mitochondrion"
        ].reset_index(drop=True)

        labels_df_filt = labels_df[
            labels_df["gene_encoding_type"] == "Mitochondrion"
        ].reset_index(drop=True)

        # Joining on the organism label data
        processed_destress_data_joined = pd.concat(
            [
                mitochondrion_df,
                labels_df_filt[
                    [
                        "organism_scientific_name",
                        "organism_group_mito",
                        # "organism_group2",
                        # "organism_group_mito",
                        # "organism_group2_mito",
                    ]
                ],
            ],
            axis=1,
        )

        # Average each principal component grouped by organism
        processed_destress_data_avg = processed_destress_data_joined.groupby(
            [
                "organism_scientific_name",
                "organism_group_mito",
            ],
            as_index=False,
        )[processed_destress_data.columns.to_list()].mean()

        # Extracting labels
        labels = processed_destress_data_avg[
            [
                "organism_scientific_name",
                # "organism_group",
                # "organism_group2",
                "organism_group_mito",
                # "organism_group2_mito",
            ]
        ]

        # Removing these labels from destress data
        processed_destress_data_avg.drop(
            [
                "organism_scientific_name",
                # "organism_group",
                # "organism_group2",
                "organism_group_mito",
                # "organism_group2_mito",
            ],
            inplace=True,
            axis=1,
        )

        # 4. Performing PCA--------------------------------------------------------------------------

        # Calculating the variance explained
        var_explained_df = pca_var_explained(
            data=processed_destress_data_avg,
            n_components=n_components,
            file_name="pca_var_explained",
            output_path=output_path_scaled,
        )

        # Performing PCA
        pca_transformed_data = perform_pca(
            data=processed_destress_data_avg,
            labels_df=labels,
            n_components=n_components,
            output_path=output_path_scaled,
            file_path="pca_transformed_data",
            components_file_path="comp_contrib",
        )

        # 5. Plotting 2d spaces---------------------------------------------------------------------

        # Setting theme for plots
        sns.set_style("whitegrid")

        # Calculating the variance explained for PC1 and PC2
        x_var_explained = var_explained_df["var_explained"][
            var_explained_df["n_components"] == 1
        ]
        y_var_explained = var_explained_df["var_explained"][
            var_explained_df["n_components"] == 2
        ]
        z_var_explained = var_explained_df["var_explained"][
            var_explained_df["n_components"] == 3
        ]

        # Formatting this for the plots
        x_var_explained_formatted = np.round(x_var_explained.iloc[0], 2) * 100
        y_var_explained_formatted = np.round(y_var_explained.iloc[0], 2) * 100
        z_var_explained_formatted = np.round(z_var_explained.iloc[0], 2) * 100

        # Scatter plot of PC1 against PC2
        plot = sns.scatterplot(
            data=pca_transformed_data,
            x="dim0",
            y="dim1",
            alpha=0.8,
            s=50,
            legend=True,
            linewidth=0.2,
            edgecolor="black",
        )
        plt.xlabel(
            "PC1 (" + str(np.int64(x_var_explained_formatted)) + "%)",
            fontsize=15,
        )
        plt.ylabel(
            "PC2 (" + str(np.int64(y_var_explained_formatted)) + "%)",
            fontsize=15,
        )
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.savefig(
            output_path_scaled + "pca_embedding_12.png",
            bbox_inches="tight",
            dpi=600,
        )
        plt.close()

        # Scatter plot of PC1 against PC3
        plot = sns.scatterplot(
            data=pca_transformed_data,
            x="dim0",
            y="dim2",
            alpha=0.8,
            s=50,
            legend=True,
            linewidth=0.2,
            edgecolor="black",
        )
        plt.xlabel(
            "PC1 (" + str(np.int64(x_var_explained_formatted)) + "%)",
            fontsize=15,
        )
        plt.ylabel(
            "PC3 (" + str(np.int64(z_var_explained_formatted)) + "%)",
            fontsize=15,
        )
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.savefig(
            output_path_scaled + "pca_embedding_13.png",
            bbox_inches="tight",
            dpi=600,
        )
        plt.close()

        # Scatter plot of PC1 against PC3
        plot = sns.scatterplot(
            data=pca_transformed_data,
            x="dim1",
            y="dim2",
            alpha=0.8,
            s=50,
            legend=True,
            linewidth=0.2,
            edgecolor="black",
        )
        plt.xlabel(
            "PC2 (" + str(np.int64(y_var_explained_formatted)) + "%)",
            fontsize=15,
        )
        plt.ylabel(
            "PC3 (" + str(np.int64(z_var_explained_formatted)) + "%)",
            fontsize=15,
        )
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.savefig(
            output_path_scaled + "pca_embedding_23.png",
            bbox_inches="tight",
            dpi=600,
        )
        plt.close()

        # Plotly scatter plot of PC1 against PC2
        fig = px.scatter(
            pca_transformed_data,
            x="dim0",
            y="dim1",
            opacity=0.9,
            hover_data=hover_data,
            labels={
                "dim0": "PC1",
                "dim1": "PC2",
            },
        )
        fig.update_traces(
            marker=dict(size=10, line=dict(width=0.8)),
            selector=dict(mode="markers"),
        )
        fig.write_html(output_path_scaled + "pca_embedding_12.html")

        plot_pca_boxplots(
            principal_components_list=["dim0", "dim1", "dim2", "dim3"],
            pca_data=pca_transformed_data,
            x="organism_group_mito",
            output_path=output_path_scaled,
            palette=palette,
        )

        # Producing the same PCA plots with the points coloured by different labels
        for label in label_dict.keys():

            cmap = palette

            hue_order = (
                pca_transformed_data.sort_values(by=label, ascending=True)[label]
                .unique()
                .tolist()
            )

            # Plotting PCA plot coloured by label
            plot_latent_space_2d(
                data=pca_transformed_data.sort_values(by=label, ascending=True),
                var_explained_data=var_explained_df,
                x="dim0",
                y="dim1",
                axes_prefix="PC",
                legend_title=label_dict[label],
                hue=label,
                hue_order=hue_order,
                # style="cluster_representative",
                alpha=0.9,
                s=140,
                palette=cmap,
                output_path=output_path_scaled,
                file_name="pca_embedding_" + label,
            )

            # Plotting PCA plot coloured by label
            plot_latent_space_2d(
                data=pca_transformed_data.sort_values(by=label, ascending=True),
                var_explained_data=var_explained_df,
                x="dim0",
                y="dim2",
                axes_prefix="PC",
                legend_title=label_dict[label],
                hue=label,
                hue_order=hue_order,
                # style="cluster_representative",
                alpha=0.9,
                s=140,
                palette=cmap,
                output_path=output_path_scaled,
                file_name="pca_embedding_" + label,
            )

            # Plotting PCA plot coloured by label
            plot_latent_space_2d(
                data=pca_transformed_data.sort_values(by=label, ascending=True),
                var_explained_data=var_explained_df,
                x="dim1",
                y="dim2",
                axes_prefix="PC",
                legend_title=label_dict[label],
                hue=label,
                hue_order=hue_order,
                # style="cluster_representative",
                alpha=0.9,
                s=140,
                palette=cmap,
                output_path=output_path_scaled,
                file_name="pca_embedding_" + label,
            )

            spectral_plot(
                pca_data=pca_transformed_data.sort_values(
                    by="organism_scientific_name", ascending=True
                ),
                group_var="organism_group_mito",
                value_var_list=dim_ids_list,
                filt_list=None,
                title="",
                legend_title="",
                output_path=output_path_scaled,
                file_name="spectral_plot",
                palette=palette,
            )

            # Plotly scatter plot of PC1 against PC2
            fig = px.scatter(
                pca_transformed_data,
                x="dim0",
                y="dim1",
                opacity=0.9,
                hover_data=hover_data,
                labels={
                    "dim0": "PC1",
                    "dim1": "PC2",
                },
                color="organism_group_mito",
                # symbol="cluster_representative",
                # color_discrete_map=palette,
            )
            fig.update_traces(
                marker=dict(size=30, line=dict(width=0.8)),
                selector=dict(mode="markers"),
            )
            fig.write_html(output_path_scaled + "pca_embedding_12.html")
