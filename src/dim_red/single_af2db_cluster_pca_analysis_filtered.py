# This script runs Principal Component Analysis (PCA)
# across the different data sets and scaling methods.

# 0. Importing packages------------------------------------------------------------
from dim_red_tools import *
import os


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

# Defining data path.
# For trimmed/aligned structures (original single-cluster workflow, A0A0G9LKG2 only):
# data_path = "data/processed_data_single_af2db_cluster/"
# Using the full processed dataset (untrimmed) as a first pass across all clusters:
data_path = "data/processed_data/"
raw_data_path = "data/raw_data/"

# Defining output data path
output_path = "analysis/pca_single_af2db_cluster_filtered/"

# Defining a list of the unipot descriptions
af2db_cluster_list_path = (
    "data/processed_data/af2/uniq_org_counts_by_structural_cluster_gt_40.csv"
)

# Defining a dictionary of labels
label_dict = {
    "organism_group": "Kingdom",
    # "organism_group2": "Domain",
}


# 2. Looping through the different data sets--------------------------------------------------

for dataset in dataset_list:

    # Defining the paths for the specified data set
    data_path_dataset = data_path + dataset + "/"
    output_path_dataset = output_path + dataset + "/"

    # 3. Building the cluster list and loading labels — done once per dataset ----------------

    # Load labels from the non-scaled path (same for all scaling methods).
    labels_df_path = data_path_dataset + "labels.csv"
    labels_df = pd.read_csv(labels_df_path)
    labels_df_for_filtering = labels_df

    # Computing unique organism counts per cluster from the pLDDT-filtered processed data,
    # so the >= 40 threshold is applied to the actual data being analysed (not pre-pLDDT counts).
    org_counts = (
        labels_df_for_filtering.groupby("cluster_representative")[
            "organism_scientific_name"
        ]
        .nunique()
        .reset_index()
        .rename(columns={"organism_scientific_name": "n_orgs"})
    )
    af2db_cluster_list_df = org_counts[org_counts["n_orgs"] >= 40].reset_index(
        drop=True
    )

    # Joining the most common UniProt description for each cluster from the labels data,
    # so we can filter on protein identity rather than hardcoding cluster IDs.
    labels_desc = labels_df_for_filtering[
        ["cluster_representative", "uniprot_description"]
    ]
    top_desc = (
        labels_desc.groupby("cluster_representative")["uniprot_description"]
        .agg(lambda x: x.value_counts().index[0])
        .reset_index()
        .rename(columns={"uniprot_description": "top_description"})
    )
    af2db_cluster_list_df = af2db_cluster_list_df.merge(
        top_desc, on="cluster_representative", how="left"
    )

    # Excluding clusters that are not well-characterised single-function proteins.
    # Two categories are excluded:
    # 1. "Uncharacterized protein" - no functional interpretation is possible so any
    #    PCA signal cannot be attributed to a specific biological process.
    # 2. Domain-containing proteins (e.g. MFS, Aminotran_5, M20_dimer, PKS_ER, Aa_trans)
    #    - these are broad structural superfamilies rather than true orthologous protein
    #    families, making it impossible to draw meaningful biological conclusions from
    #    the PCA.
    exclude_terms = ["Uncharacterized protein", "domain-containing protein"]
    exclude_mask = af2db_cluster_list_df["top_description"].str.contains(
        "|".join(exclude_terms), na=False
    )
    af2db_cluster_list_df = af2db_cluster_list_df[~exclude_mask].reset_index(drop=True)

    af2db_cluster_list = af2db_cluster_list_df["cluster_representative"].to_list()

    print(f"Running {len(af2db_cluster_list)} clusters after filtering:")
    for c in af2db_cluster_list:
        desc = af2db_cluster_list_df.loc[
            af2db_cluster_list_df["cluster_representative"] == c, "top_description"
        ].values[0]
        print(f"  {c}: {desc}")

    # To run on all 22 clusters including uncharacterised and domain-containing proteins:
    # af2db_cluster_list = pd.read_csv(af2db_cluster_list_path)["cluster_representative"].to_list()
    # To run on a single cluster for testing:
    # af2db_cluster_list = ["A0A0G9LKG2"]

    for scaling_method in scaling_method_list:

        # Defining paths for the data and scaling method used
        data_path_scaled = data_path_dataset + scaling_method + "/"
        output_path_scaled = output_path_dataset + scaling_method + "/"

        # Defining the path for processed AF2 DE-STRESS data
        processed_destress_data_path = (
            data_path_scaled + "processed_destress_data_scaled.csv"
        )

        # 3. Reading in data------------------------------------------------------------------------

        # Defining the path for the processed AF2 DE-STRESS data
        processed_destress_data = pd.read_csv(processed_destress_data_path)

        # labels_df is loaded once outside the scaling method loop and reused here.

        # processed_destress_data = processed_destress_data[
        #     ~labels_df["design_name"].isin(
        #         ["AF-Q7HP05-F1-model_v4", "AF-A0A5K1K958-F1-model_v4"]
        #     )
        # ].reset_index(drop=True)

        # labels_df = labels_df[
        #     ~labels_df["design_name"].isin(
        #         ["AF-Q7HP05-F1-model_v4", "AF-A0A5K1K958-F1-model_v4"]
        #     )
        # ].reset_index(drop=True)

        for af2db_cluster in af2db_cluster_list:

            # Filtering destress data
            processed_destress_data_filt = processed_destress_data[
                labels_df["cluster_representative"] == af2db_cluster
            ].reset_index(drop=True)

            # Filtering labels data
            labels_df_filt = labels_df[
                labels_df["cluster_representative"] == af2db_cluster
            ].reset_index(drop=True)

            # Creating a subfolder for this uniprot description in the output directory
            create_subfolder(output_path_scaled, af2db_cluster)

            # Creating new output directory
            output_path_af2db_cluster = output_path_scaled + af2db_cluster + "/"

            # 4. Performing PCA--------------------------------------------------------------------------

            # Calculating the variance explained
            var_explained_df = pca_var_explained(
                data=processed_destress_data_filt,
                n_components=n_components,
                file_name="pca_var_explained",
                output_path=output_path_af2db_cluster,
            )

            # Performing PCA
            pca_transformed_data = perform_pca(
                data=processed_destress_data_filt,
                labels_df=labels_df_filt,
                n_components=n_components,
                output_path=output_path_af2db_cluster,
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
                output_path_af2db_cluster + "pca_embedding_12.png",
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
                output_path_af2db_cluster + "pca_embedding_13.png",
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
                output_path_af2db_cluster + "pca_embedding_23.png",
                bbox_inches="tight",
                dpi=600,
            )
            plt.close()

            plot_pca_boxplots(
                principal_components_list=["dim0", "dim1", "dim2", "dim3"],
                pca_data=pca_transformed_data,
                x="organism_group",
                output_path=output_path_af2db_cluster,
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

                # hue_order = [
                #     "Animal",
                #     "Archaea",
                #     "Bacteria",
                #     "Fungi",
                #     "Plant",
                #     "Protozoan",
                # ]

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
                    style=None,
                    alpha=0.9,
                    s=140,
                    palette=cmap,
                    output_path=output_path_af2db_cluster,
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
                    style=None,
                    alpha=0.9,
                    s=140,
                    palette=cmap,
                    output_path=output_path_af2db_cluster,
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
                    style=None,
                    alpha=0.9,
                    s=140,
                    palette=cmap,
                    output_path=output_path_af2db_cluster,
                    file_name="pca_embedding_" + label,
                )

                # spectral_plot(
                #     pca_data=pca_transformed_data.sort_values(
                #         by="organism_scientific_name", ascending=True
                #     ),
                #     group_var="organism_group_mito",
                #     value_var_list=dim_ids_list,
                #     filt_list=None,
                #     title=af2db_cluster,
                #     legend_title="",
                #     output_path=output_path_af2db_cluster,
                #     file_name="spectral_plot",
                #     palette=palette,
                # )

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
                    # symbol="cluster_representative",
                    # color_discrete_map=palette,
                )
                fig.update_traces(
                    marker=dict(size=30, line=dict(width=0.8)),
                    selector=dict(mode="markers"),
                )
                fig.write_html(output_path_af2db_cluster + "pca_embedding_12.html")
