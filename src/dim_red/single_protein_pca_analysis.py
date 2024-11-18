# This script runs Principal Component Analysis (PCA)
# across the different data sets and scaling methods.

# 0. Importing packages------------------------------------------------------------
from dim_red_tools import *
import os


def dataframe_to_fasta(df, filename="output.fasta"):
    """
    Converts a DataFrame with columns 'design_name' and 'full_sequence' to a .fasta format file.

    Parameters:
    - df: pandas DataFrame containing the data.
    - filename: File name where the FASTA data will be written.
    """
    with open(filename, "w") as file:
        for _, row in df.iterrows():
            fasta_format = f">{row['design_name']}\n{row['full_sequence']}\n"
            file.write(fasta_format)
    print(f"FASTA file created: {filename}")


def create_subfolder(directory, subfolder):
    path = os.path.join(directory, subfolder)
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Directory '{path}' created")
    else:
        print(f"Directory '{path}' already exists")


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
hover_data = ["design_name", "dim0", "dim1", "organism_scientific_name"]

# Creating a color palette
palette = sns.color_palette(
    ["#0173b2", "#d55e00", "#029e73", "#cc78bc", "#808080", "#f0e442"], 6
)

# palette = sns.color_palette("tab10")

# Defining data path
data_path = "data/processed_data/"
raw_data_path = "data/raw_data/"

# Defining output data path
output_path = "analysis/pca_single_proteins/"

# Defining a list of the unipot descriptions
uniprot_desc_list_path = (
    "analysis/data_exploration/uniprot/uniprot_desc_org_count_gt_40.csv"
)

# Defining a dictionary of labels
label_dict = {
    "organism_group_mito": "Organism Group 1",
    "organism_group2_mito": "Organism Group 2",
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

        # # Filtering destress data
        # processed_destress_data = processed_destress_data[
        #     (labels_df["Mean_PLDDT"] >= 80)
        #     & ~(
        #         labels_df["design_name"].isin(
        #             ["AF-Q5F599-F1-model_v4", "AF-K7L111-F1-model_v4"]
        #         )
        #     )
        # ].reset_index(drop=True)

        # # Filtering labels data
        # labels_df = labels_df[
        #     (labels_df["Mean_PLDDT"] >= 80)
        #     & ~labels_df["design_name"].isin(
        #         ["AF-A0A0R0KJP6-F1-model_v4", "AF-A0A3Q0KVD7-F1-model_v4"]
        #     )
        # ].reset_index(drop=True)

        # Filtering destress data
        processed_destress_data = processed_destress_data[
            (labels_df["Mean_PLDDT"] >= 80)
            # & ~(
            #     labels_df["design_name"].isin(
            #         ["AF-A0A0R0KJP6-F1-model_v4", "AF-A0A3Q0KVD7-F1-model_v4"]
            #     )
            # )
        ].reset_index(drop=True)

        # Filtering labels data
        labels_df = labels_df[
            (labels_df["Mean_PLDDT"] >= 80)
            # & ~labels_df["design_name"].isin(
            #     ["AF-Q5F599-F1-model_v4", "AF-K7L111-F1-model_v4"]
            # )
        ].reset_index(drop=True)

        # Filtering for gene_encoding_type 'Mitochondrion'
        mitochondrion_df = labels_df[labels_df["gene_encoding_type"] == "Mitochondrion"]

        # Grouping by 'uniprot_description' and counting unique 'organism_scientific_name'
        grouped_descriptions = mitochondrion_df.groupby("uniprot_description")[
            "organism_scientific_name"
        ].nunique()

        # Filtering descriptions with 10 or more unique organisms
        descriptions_with_10plus_organisms = grouped_descriptions[
            grouped_descriptions >= 10
        ]

        # Listing descriptions
        uniprot_desc_list = descriptions_with_10plus_organisms.index.tolist()

        print(uniprot_desc_list)

        # # Reading in uniprot description list
        # uniprot_desc_list = pd.read_csv(uniprot_desc_list_path)
        # uniprot_desc_list = uniprot_desc_list["uniprot_description"].to_list()

        uniprot_desc_list = ["NADH-ubiquinone oxidoreductase chain 4"]
        # uniprot_desc_list = ["Cytochrome b"]
        # uniprot_desc_list = ["ATP synthase subunit a"]
        # U5E5N1

        # processed_destress_data = processed_destress_data[
        #     labels_df["uniprot_description"] == "tRNA (guanine-N(7)-)-methyltransferase"
        # ].reset_index(drop=True)

        # labels_df = labels_df[
        #     labels_df["uniprot_description"] == "tRNA (guanine-N(7)-)-methyltransferase"
        # ].reset_index(drop=True)

        # labels_df.to_csv("trna_methyl_transferase.csv", index=False)

        # # Convert to FASTA format and print
        # dataframe_to_fasta(labels_df[["design_name", "full_sequence"]])

        for uniprot_description in uniprot_desc_list:
            if uniprot_description != "Uncharacterized protein":

                # Filtering destress data
                processed_destress_data_filt = processed_destress_data[
                    (labels_df["uniprot_description"] == uniprot_description)
                    # & (labels_df["Mean_PLDDT"] >= 80)
                ].reset_index(drop=True)

                # Filtering labels data
                labels_df_filt = labels_df[
                    (labels_df["uniprot_description"] == uniprot_description)
                    # & (labels_df["Mean_PLDDT"] >= 80)
                ].reset_index(drop=True)

                # print(
                #     labels_df_filt[["design_name", "uniprot_description", "Mean_PLDDT"]]
                # )

                # labels_df_filt["mito_flag"] = labels_df_filt[
                #     "subcellular_location"
                # ].str.contains("Mitochondrion")

                # labels_df_filt["organism_group"] = np.where(
                #     labels_df_filt["mito_flag"] == 1,
                #     labels_df_filt["organism_group"]
                #     + "- Mitochondrion",  # value when mito_flag is 1
                #     labels_df_filt["organism_group"],  # keep original value otherwise
                # )

                # processed_destress_data_filt = processed_destress_data[
                #     labels_df["cluster_representative"] == "A0A4W4EL89"
                # ].reset_index(drop=True)

                # labels_df_filt = labels_df[
                #     labels_df["cluster_representative"] == "A0A4W4EL89"
                # ].reset_index(drop=True)

                # "ATP synthase subunit a"

                # Creating a subfolder for this uniprot description in the output directory
                create_subfolder(output_path_scaled, uniprot_description)

                # Creating new output directory
                output_path_uniprot_desc = (
                    output_path_scaled + uniprot_description + "/"
                )

                # 4. Performing PCA--------------------------------------------------------------------------

                # Calculating the variance explained
                var_explained_df = pca_var_explained(
                    data=processed_destress_data_filt,
                    n_components=n_components,
                    file_name="pca_var_explained",
                    output_path=output_path_uniprot_desc,
                )

                # Performing PCA
                pca_transformed_data = perform_pca(
                    data=processed_destress_data_filt,
                    labels_df=labels_df_filt,
                    n_components=n_components,
                    output_path=output_path_uniprot_desc,
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
                    output_path_uniprot_desc + "pca_embedding_12.png",
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
                    output_path_uniprot_desc + "pca_embedding_13.png",
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
                    output_path_uniprot_desc + "pca_embedding_23.png",
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
                fig.write_html(output_path_uniprot_desc + "pca_embedding_12.html")

                plot_pca_boxplots(
                    principal_components_list=["dim0", "dim1", "dim2", "dim3"],
                    pca_data=pca_transformed_data,
                    x="organism_group_mito",
                    output_path=output_path_uniprot_desc,
                    palette=palette,
                )

                # Producing the same PCA plots with the points coloured by different labels
                for label in label_dict.keys():

                    cmap = palette

                    hue_order = (
                        pca_transformed_data.sort_values(by=label, ascending=True)[
                            label
                        ]
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
                        output_path=output_path_uniprot_desc,
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
                        output_path=output_path_uniprot_desc,
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
                        output_path=output_path_uniprot_desc,
                        file_name="pca_embedding_" + label,
                    )

                    spectral_plot(
                        pca_data=pca_transformed_data.sort_values(
                            by="organism_scientific_name", ascending=True
                        ),
                        group_var="organism_group",
                        value_var_list=dim_ids_list,
                        filt_list=None,
                        title=uniprot_description,
                        legend_title="",
                        output_path=output_path_uniprot_desc,
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
                        symbol="cluster_representative",
                        # color_discrete_map=palette,
                    )
                    fig.update_traces(
                        marker=dict(size=30, line=dict(width=0.8)),
                        selector=dict(mode="markers"),
                    )
                    fig.write_html(output_path_uniprot_desc + "pca_embedding_12.html")
