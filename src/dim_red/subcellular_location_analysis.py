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

# Defining a dictionary of labels
label_dict = {
    "organism_group_mito": "Organism Group 1",
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
    ].reset_index(drop=True)

    # Filtering labels data
    labels_df = labels_df[(labels_df["Mean_PLDDT"] >= 70)].reset_index(drop=True)

    # Outputting unique values for subcellular location acrossd the data set
    uniq_subceullar_location = labels_df["subcellular_location"].value_counts()
    uniq_subceullar_location.to_csv(
        output_path_scaled + "subcellular_location_unique_value_count.csv"
    )

    # Outputting all the singleton subcellular locations
    labels_df["subcellular_location"] = labels_df["subcellular_location"].astype(str)

    # Generate the mask, where NaN in the original data leads to False
    mask = labels_df["subcellular_location"].str.contains(":").fillna(False)

    # Apply the NOT operator
    final_mask = ~mask

    # Now you can use final_mask for further processing or filtering
    labels_df_filt = labels_df[final_mask]

    # Outputting unique values for subcellular location acrossd the data set
    uniq_subceullar_location_single = labels_df_filt[
        "subcellular_location"
    ].value_counts()
    uniq_subceullar_location_single.to_csv(
        output_path_scaled + "subcellular_location_unique_value_count_single.csv"
    )

    top_subceullar_locations = [
        "Membrane",
        "Nucleus",
        "Cytoplasm",
        "Cell membrane",
        "Secreted",
        "Mitochondrion",
        "Cell inner membrane",
        "Endoplasmic reticulum membrane",
        "Mitochondrion inner membrane",
        "Nucleus, nucleolus",
        "Plastid, chloroplast",
        "Cytoplasm, cytoskeleton",
        "Golgi apparatus membrane",
    ]

    uniq_subceullar_location_single_filt = uniq_subceullar_location_single.loc[
        top_subceullar_locations
    ]

    # Plotting the data
    plt.figure(figsize=(10, 8))  # Adjust the size of the figure as needed
    uniq_subceullar_location_single_filt.plot(
        kind="bar", color="blue", x="subcellular_location", y="count"
    )  # You can customize the color
    plt.title(f"Protein count per subcellular location")
    plt.ylabel("Protein Count")
    plt.xlabel("Subcellular Location")
    plt.xticks(rotation=45, ha="right")  # Rotate labels for better readability
    plt.tight_layout()  # Adjust subplots to fit into figure area.

    # Save the figure
    plot_filename = "subceullar_location_protein_count.png"
    plt.savefig(os.path.join(output_path_scaled, plot_filename))
    plt.close()  # Close the plot to free up memory

    # Filtering all the labels by the top subcellular locations
    labels_df_top_subcellular_locations = labels_df[
        labels_df["subcellular_location"].isin(top_subceullar_locations)
    ].reset_index(drop=True)

    labels_df_top_subcellular_locations = labels_df_top_subcellular_locations[
        "organism_scientific_name"
    ].value_counts()

    # Plotting the data
    plt.figure(figsize=(10, 8))  # Adjust the size of the figure as needed
    labels_df_top_subcellular_locations.plot(
        kind="bar", color="blue", x="organism_scientific_name", y="count"
    )  # You can customize the color
    plt.title(f"Protein count per organism")
    plt.ylabel("Protein Count")
    plt.xlabel("Organism")
    plt.xticks(rotation=45, ha="right")  # Rotate labels for better readability
    plt.tight_layout()  # Adjust subplots to fit into figure area.

    # Save the figure
    plot_filename = "organism_protein_count_top_subceullular_locations.png"
    plt.savefig(os.path.join(output_path_scaled, plot_filename))
    plt.close()  # Close the plot to free up memory
