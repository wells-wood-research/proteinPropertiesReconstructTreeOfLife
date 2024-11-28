# 0. Importing packages----------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn import metrics
from clustering_tools import *

# 1. Defining variables-----------------------------------------------------------------

# Defining the scaling methods list
scaling_method_list = ["standard", "robust", "minmax"]
# scaling_method_list = ["standard"]

# Creating a color palette
palette = sns.color_palette(
    ["#0173b2", "#d55e00", "#029e73", "#cc78bc", "#808080", "#f0e442"], 6
)

# Defining data path
data_path = "data/processed_data/af2/"
raw_data_path = "data/raw_data/"

# Defining output data path
output_path = "analysis/kmeans_avg_by_org/af2/"

# Defining a dictionary of labels
label_dict = {
    "organism_group": "Organism Group 1",
    "organism_group2": "Organism Group 2",
}

# Defining the number of clusters
n_clusters_list = range(2, 21, 1)

# Defining the number of kmeans initialisations
n_inits = 100

# Creating a data frame to gather the clustering results
clustering_results_master = pd.DataFrame(
    columns=[
        "model",
        "dataset",
        "scaler",
        "kmeans_init",
        "n_clusters",
        "weighted_ssd",
        "adj_rand_score",
    ]
)

# Grouping var
group_var = "organism_group2"

# 2. Looping through the different scaling methods--------------------------------------------------

for scaling_method in scaling_method_list:

    # Defining paths for the data and scaling method used
    data_path_scaled = data_path + scaling_method + "/"
    output_path_scaled = output_path + scaling_method + "/"

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

    # 4. Running different initialisations of k means--------------------------------------------

    for n_clusters in n_clusters_list:
        for init in range(0, n_inits, 1):
            # Generating a random integer
            rand_int = np.random.randint(0, high=100000, size=1)[0]

            # Setting up kmeans
            model = KMeans(
                n_clusters=n_clusters,
                n_init="auto",
                random_state=rand_int,
            )
            model_fit = model.fit(processed_destress_data_avg)

            # Extracting the labels
            predicted_labels = model_fit.labels_

            # Extracting the sum of squared distances of samples to
            # their closest cluster centre
            weighted_ssd = model.inertia_

            # Calculating the adjusted rand score against
            # the organism labels
            adj_rand_score = metrics.adjusted_rand_score(
                eval(group_var + "_labels"),
                predicted_labels,
            )

            # Creating a row data frame
            clustering_results = pd.DataFrame(
                {
                    "model": "kmeans",
                    "dataset": "af2",
                    "scaler": scaling_method,
                    "kmeans_init": init,
                    "n_clusters": n_clusters,
                    "weighted_ssd": weighted_ssd,
                    "adj_rand_score": adj_rand_score,
                },
                index=[0],
            )

            # Adding the hyper parameters to the data set
            clustering_results_master = pd.concat(
                [clustering_results_master, clustering_results],
                axis=0,
                ignore_index=True,
            )

clustering_results_master.to_csv(
    output_path + "kmeans_results_master_destress" + "_" + group_var + ".csv",
    index=False,
)

# Plotting the average weighted_ssd and adj_rand_score by scaler and number of clusters
# for k means and adj_rand_score by scaler and number of clusters for hierarchical clustering
for scaling_method in scaling_method_list:

    # Defining paths for the data and scaling method used
    output_path_scaled = output_path + scaling_method + "/"

    clustering_results_master_scaler = clustering_results_master[
        clustering_results_master["scaler"] == scaling_method
    ].reset_index(drop=True)

    plt.figure(figsize=(6, 5))
    sns.set_style("whitegrid")

    adj_rand_ind_plot(
        data=clustering_results_master_scaler,
        title="",
        file_name="kmeans_eval_" + scaling_method + "_destress_" + group_var,
        # hue="group_var",
        output_path=output_path_scaled,
    )
