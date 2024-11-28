# This script runs Principal Component Analysis (PCA)
# across the different data sets and scaling methods.

# 0. Importing packages------------------------------------------------------------
from dim_red_tools import *
import os

# 1. Defining variables------------------------------------------------------------

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
output_path = "analysis/pca_subcellular_location/"
# output_path = "analysis/pca_subcellular_location_euk_cyto/"

# Defining the list of subcellular locations
subcellular_locations = [
    "Membrane",
    "Nucleus",
    "Cytoplasm",
]

# Defining a dictionary of labels
label_dict = {
    "Subcellular Location": "Subcellular Location & Organism Kingdom",
    # "Kingdom": "Kingdom",
    # "organism_group2_mito": "Organism Group 2",
}

for scaling_method in scaling_method_list:

    # Defining paths for the data and scaling method used
    data_path_scaled = data_path + "af2/" + scaling_method + "/"
    output_path_scaled = output_path + "af2/" + scaling_method + "/"

    # Defining the path for processed AF2 DE-STRESS data
    processed_destress_data_path = (
        data_path_scaled + "processed_destress_data_scaled.csv"
    )

    # Defining file paths for labels
    labels_df_path = data_path + "af2/" + "labels.csv"

    # 3. Reading in data------------------------------------------------------------------------

    # Defining the path for the processed AF2 DE-STRESS data
    processed_destress_data = pd.read_csv(processed_destress_data_path)

    # Reading in labels
    labels_df = pd.read_csv(labels_df_path)

    # Filtering destress data
    processed_destress_data = processed_destress_data[
        (labels_df["Mean_PLDDT"] >= 70)
        & (labels_df["subcellular_location"].isin(subcellular_locations))
        & (labels_df["organism_group2"] == "Eukaryotes")
    ].reset_index(drop=True)

    # Filtering labels data
    labels_df = labels_df[
        (labels_df["Mean_PLDDT"] >= 70)
        & (labels_df["subcellular_location"].isin(subcellular_locations))
        & (labels_df["organism_group2"] == "Eukaryotes")
    ].reset_index(drop=True)

    # Joining on the organism label data
    processed_destress_data_joined = pd.concat(
        [
            processed_destress_data,
            labels_df[
                [
                    "organism_scientific_name",
                    "organism_group",
                    "subcellular_location",
                ]
            ],
        ],
        axis=1,
    )

    # Average each principal component grouped by organism
    processed_destress_data_avg = processed_destress_data_joined.groupby(
        [
            "organism_scientific_name",
            "organism_group",
            "subcellular_location",
        ],
        as_index=False,
    )[processed_destress_data.columns.to_list()].mean()

    # Extracting labels
    labels = processed_destress_data_avg[
        [
            "organism_scientific_name",
            "organism_group",
            "subcellular_location",
        ]
    ]

    # Removing these labels from destress data
    processed_destress_data_avg.drop(
        [
            "organism_scientific_name",
            "organism_group",
            "subcellular_location",
        ],
        inplace=True,
        axis=1,
    )

    # Rename columns
    labels.rename(
        columns={
            "organism_group2": "Domain",
            "organism_group": "Kingdom",
            "subcellular_location": "Subcellular Location",
        },
        inplace=True,
    )

    print(labels)

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
        x="Kingdom",
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
            style="Kingdom",
            # style=None,
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
            style="Kingdom",
            # style=None,
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
            style="Kingdom",
            # style=None,
            alpha=0.9,
            s=140,
            palette=cmap,
            output_path=output_path_scaled,
            file_name="pca_embedding_" + label,
        )

        # # Plotly scatter plot of PC1 against PC2
        # fig = px.scatter(
        #     pca_transformed_data,
        #     x="dim0",
        #     y="dim1",
        #     opacity=0.9,
        #     hover_data=hover_data,
        #     labels={
        #         "dim0": "PC1",
        #         "dim1": "PC2",
        #     },
        #     color="Kingdom",
        #     # color_discrete_map=palette,
        # )
        # fig.update_traces(
        #     marker=dict(size=30, line=dict(width=0.8)),
        #     selector=dict(mode="markers"),
        # )
        # fig.write_html(output_path_scaled + "pca_embedding_12.html")
