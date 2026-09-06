# K-means baseline using amino acid composition (20 standard AAs).
# Reads scaled AA composition data produced by data_prep_aa_composition.py,
# averages per organism, and evaluates k = 2-20 with 100 initialisations each.
# Outputs clustering results CSV and adjusted Rand index plots to
# analysis/kmeans_avg_by_org_aa_composition/af2/.

# 0. Importing packages---------------------------------------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn import metrics
from clustering_tools import *

# 1. Defining variables---------------------------------------------------------

scaling_method_list = ["standard", "robust", "minmax"]

data_path = "data/processed_data_aa_composition/af2/"

output_path = "analysis/kmeans_avg_by_org_aa_composition/af2/"

n_clusters_list = range(2, 21, 1)

n_inits = 100

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

# 2. Looping through the different scaling methods------------------------------

for scaling_method in scaling_method_list:

    data_path_scaled = data_path + scaling_method + "/"
    output_path_scaled = output_path + scaling_method + "/"

    scaled_data_path = data_path_scaled + "aa_composition_scaled_nonredundant.csv"
    labels_df_path = data_path + "labels_nonredundant.csv"

    # 3. Reading in data--------------------------------------------------------

    scaled_data = pd.read_csv(scaled_data_path)
    labels_df = pd.read_csv(labels_df_path)

    data_joined = pd.concat(
        [
            scaled_data,
            labels_df[
                [
                    "organism_scientific_name",
                    "organism_group",
                    "organism_group2",
                ]
            ],
        ],
        axis=1,
    )

    feature_cols = scaled_data.columns.to_list()

    data_avg = data_joined.groupby(
        [
            "organism_scientific_name",
            "organism_group",
            "organism_group2",
        ],
        as_index=False,
    )[feature_cols].mean()

    organism_group_labels = data_avg["organism_group"].to_list()
    organism_labels = data_avg["organism_scientific_name"].to_list()

    data_avg.drop(
        [
            "organism_scientific_name",
            "organism_group",
            "organism_group2",
        ],
        inplace=True,
        axis=1,
    )

    # 4. Running k-means evaluations--------------------------------------------

    for n_clusters in n_clusters_list:
        for init in range(0, n_inits, 1):

            rand_int = np.random.randint(0, high=100000, size=1)[0]

            model = KMeans(
                n_clusters=n_clusters,
                n_init="auto",
                random_state=rand_int,
            )
            model_fit = model.fit(data_avg)

            predicted_labels = model_fit.labels_

            weighted_ssd = model.inertia_

            adj_rand_score = metrics.adjusted_rand_score(
                organism_group_labels,
                predicted_labels,
            )

            clustering_results = pd.DataFrame(
                {
                    "model": "kmeans",
                    "dataset": "af2_aa_composition",
                    "scaler": scaling_method,
                    "kmeans_init": init,
                    "n_clusters": n_clusters,
                    "weighted_ssd": weighted_ssd,
                    "adj_rand_score": adj_rand_score,
                },
                index=[0],
            )

            clustering_results_master = pd.concat(
                [clustering_results_master, clustering_results],
                axis=0,
                ignore_index=True,
            )

clustering_results_master.to_csv(
    output_path + "kmeans_results_master_aa_composition_organism_group.csv",
    index=False,
)

# 5. Plotting adjusted Rand index per scaler------------------------------------

for scaling_method in scaling_method_list:

    output_path_scaled = output_path + scaling_method + "/"

    clustering_results_master_scaler = clustering_results_master[
        clustering_results_master["scaler"] == scaling_method
    ].reset_index(drop=True)

    plt.figure(figsize=(6, 5))
    sns.set_style("whitegrid")

    adj_rand_ind_plot(
        data=clustering_results_master_scaler,
        title="",
        file_name="kmeans_eval_" + scaling_method + "_aa_composition_organism_group",
        output_path=output_path_scaled,
    )
