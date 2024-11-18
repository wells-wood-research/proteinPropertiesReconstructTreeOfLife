# This script runs Principal Component Analysis (PCA)
# across the different data sets and scaling methods.

# 0. Importing packages------------------------------------------------------------
from dim_red_tools import *

# 1. Defining variables------------------------------------------------------------

# Defining the data set list
dataset_list = ["pdb"]

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
hover_data = ["design_name", "dim0", "dim1"]

# Creating a color palette
palette = sns.color_palette(
    ["#0173b2", "#d55e00", "#029e73", "#cc78bc", "#de8f05", "#6a3d9a"], 4
)

# Defining data path
data_path = "data/processed_data/"

# Defining output data path
output_path = "analysis/pca_all_af2_models/"

# Defining a dictionary of labels
label_dict = {
    "isoelectric_point": "Isoelectric Point",
    "packing_density": "Packing Density",
    "aggrescan3d_avg_value": "Aggrescan3D Average Value",
    "dssp_bin": "Secondary Structure",
}

# Defining lsits of organisms for the spectral plots
organism_animal_list = [
    "Danio rerio",
    "Mus musculus",
    "Rattus norvegicus",
    "Homo sapiens",
]

organism_fungi_list = [
    "Candida albicans",
    "Saccharomyces cerevisiae",
    "Ajellomyces capsulatus",
    "Schizosaccharomyces pombe",
]
organism_bacteria_list = [
    "Escherichia coli",
    "Mycobacterium leprae",
    "Mycobacterium tuberculosis",
    "Salmonella typhimurium",
]
organism_plant_list = [
    "Arabidopsis thaliana",
    "Glycine max",
    "Oryza sativa",
    "Zea mays",
]
organism_protozoan_list = [
    "Plasmodium falciparum",
    "Dictyostelium discoideum",
    "Leishmania infantum",
    "Trypanosoma brucei",
    "Trypanosoma cruzi",
]


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

        # 4. Performing PCA--------------------------------------------------------------------------

        # Calculating the variance explained
        var_explained_df = pca_var_explained(
            data=processed_destress_data,
            n_components=n_components,
            file_name="pca_var_explained",
            output_path=output_path_scaled,
        )

        # Performing PCA
        pca_transformed_data = perform_pca(
            data=processed_destress_data,
            labels_df=labels_df,
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

        # Formatting this for the plots
        x_var_explained_formatted = np.round(x_var_explained.iloc[0], 2) * 100
        y_var_explained_formatted = np.round(y_var_explained.iloc[0], 2) * 100

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
            "PC1 (" + str(np.int64(x_var_explained_formatted)) + "%)", fontsize=15
        )
        plt.ylabel(
            "PC2 (" + str(np.int64(y_var_explained_formatted)) + "%)", fontsize=15
        )
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.savefig(
            output_path_scaled + "pca_embedding_12.png",
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

        # Producing the same PCA plots with the points coloured by different labels
        for label in label_dict.keys():
            if label in [
                "isoelectric_point",
                "packing_density",
                "hydrophobic_fitness",
                "aggrescan3d_avg_value",
            ]:
                cmap = sns.color_palette("viridis", as_cmap=True)

            else:
                cmap = palette

            if label == "dssp_bin":
                hue_order = ["Alpha Helix", "Beta Strand", "Loop", "Mixed"]
            else:
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
                # style=var,
                alpha=0.9,
                s=70,
                palette=cmap,
                output_path=output_path_scaled,
                file_name="pca_embedding_" + label,
            )

            if dataset == "af2":

                # Computing spectral plots for different organisms
                spectral_plot(
                    pca_data=pca_transformed_data.sort_values(
                        by="organism_scientific_name", ascending=True
                    ),
                    group_var="organism_scientific_name",
                    value_var_list=dim_ids_list,
                    filt_list=organism_plant_list,
                    title="Plant",
                    legend_title="",
                    output_path=output_path_scaled,
                    file_name="spectral_plot_plant",
                    palette=palette,
                )

                spectral_plot(
                    pca_data=pca_transformed_data.sort_values(
                        by="organism_scientific_name", ascending=True
                    ),
                    group_var="organism_scientific_name",
                    value_var_list=dim_ids_list,
                    filt_list=organism_bacteria_list,
                    title="Bacteria",
                    legend_title="",
                    output_path=output_path_scaled,
                    file_name="spectral_plot_bacteria",
                    palette=palette,
                )

                spectral_plot(
                    pca_data=pca_transformed_data.sort_values(
                        by="organism_scientific_name", ascending=True
                    ),
                    group_var="organism_scientific_name",
                    value_var_list=dim_ids_list,
                    filt_list=organism_animal_list,
                    title="Animal",
                    legend_title="",
                    output_path=output_path_scaled,
                    file_name="spectral_plot_animal",
                    palette=palette,
                )

                spectral_plot(
                    pca_data=pca_transformed_data.sort_values(
                        by="organism_scientific_name", ascending=True
                    ),
                    group_var="organism_scientific_name",
                    value_var_list=dim_ids_list,
                    filt_list=organism_fungi_list,
                    title="Funghi",
                    legend_title="",
                    output_path=output_path_scaled,
                    file_name="spectral_plot_funghi",
                    palette=palette,
                )

                spectral_plot(
                    pca_data=pca_transformed_data.sort_values(
                        by="organism_scientific_name", ascending=True
                    ),
                    group_var="organism_scientific_name",
                    value_var_list=dim_ids_list,
                    filt_list=organism_protozoan_list,
                    title="Protozoan",
                    legend_title="",
                    output_path=output_path_scaled,
                    file_name="spectral_plot_protozoan",
                    palette=palette,
                )
