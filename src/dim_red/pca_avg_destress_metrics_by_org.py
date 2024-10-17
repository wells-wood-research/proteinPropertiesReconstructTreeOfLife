# This script performs PCA on the average DE-STRESS metrics by organism

# 0. Importing packages and custom functions--------------------------------
from dim_red_tools import *

# 1. Defining variables-----------------------------------------------------

# Defining the data set list
dataset_list = ["af2"]

# Defining the scaling methods list
scaling_method_list = ["standard", "robust", "minmax"]
# scaling_method_list = ["standard"]

# Defining number of principal components
n_components = 7

# Defining list of dim ids
dim_ids_list = []
for i in range(0, n_components):
    dim_ids_list.append("dim" + str(i))


# Defining hiver data for plotly
hover_data = ["dim0", "dim1", "organism_scientific_name"]

# Creating a color palette
palette = sns.color_palette(
    ["#0173b2", "#d55e00", "#029e73", "#cc78bc", "#808080", "#f0e442"], 6
)

# palette = sns.color_palette("tab10")

# Defining data path
data_path = "data/processed_data/"
raw_data_path = "data/raw_data/"

# Defining output data path
output_path = "analysis/pca_avg_by_org/"


# Defining a dictionary of labels
label_dict = {
    "organism_group": "Organism Group 1",
    "organism_group2": "Organism Group 2",
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

        # Average each principal component grouped by organism
        processed_destress_data_avg = processed_destress_data_joined.groupby(
            ["organism_scientific_name", "organism_group", "organism_group2"],
            as_index=False,
        )[processed_destress_data.columns.to_list()].median()

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

        # 5. Plotting PCA spaces---------------------------------------------------------------------

        plot_pca_boxplots(
            principal_components_list=["dim0", "dim1", "dim2", "dim3"],
            pca_data=pca_transformed_data,
            x="organism_group",
            output_path=output_path_scaled,
            palette=palette,
        )

        # Producing the same PCA plots with the points coloured by different labels
        for label in label_dict.keys():

            cmap = palette

            hue_order = (
                pca_transformed_data.sort_values(by=label, ascending=False)[label]
                .unique()
                .tolist()
            )

            # Plotting PCA plot coloured by label
            plot_latent_space_2d(
                data=pca_transformed_data.sort_values(by=label, ascending=False),
                var_explained_data=var_explained_df,
                x="dim0",
                y="dim1",
                axes_prefix="PC",
                legend_title=label_dict[label],
                hue=label,
                hue_order=hue_order,
                alpha=0.9,
                s=140,
                palette=cmap,
                output_path=output_path_scaled,
                file_name="pca_embedding_" + label,
            )

            # Plotting PCA plot coloured by label
            plot_latent_space_2d(
                data=pca_transformed_data.sort_values(by=label, ascending=False),
                var_explained_data=var_explained_df,
                x="dim0",
                y="dim2",
                axes_prefix="PC",
                legend_title=label_dict[label],
                hue=label,
                hue_order=hue_order,
                alpha=0.9,
                s=140,
                palette=cmap,
                output_path=output_path_scaled,
                file_name="pca_embedding_" + label,
            )

            # Plotting PCA plot coloured by label
            plot_latent_space_2d(
                data=pca_transformed_data.sort_values(by=label, ascending=False),
                var_explained_data=var_explained_df,
                x="dim1",
                y="dim2",
                axes_prefix="PC",
                legend_title=label_dict[label],
                hue=label,
                hue_order=hue_order,
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
                group_var="organism_group",
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
                color="organism_group",
                # color_discrete_map=palette,
            )
            fig.update_traces(
                marker=dict(size=30, line=dict(width=0.8)),
                selector=dict(mode="markers"),
            )
            fig.write_html(output_path_scaled + "pca_embedding_12.html")
