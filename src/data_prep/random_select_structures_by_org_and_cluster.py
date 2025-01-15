import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def select_random_by_organism_and_cluster_optimized(df, seed):
    """
    Optimized function to randomly select one row for each combination of 'organism_scientific_name'
    and 'cluster_representative' from a given DataFrame, using a random seed for reproducibility.

    Parameters:
    df (pandas.DataFrame): DataFrame with at least columns 'organism_scientific_name' and 'cluster_representative'.
    seed (int): Random seed for reproducibility.

    Returns:
    pandas.DataFrame: A DataFrame with randomly selected rows for each unique combination of
                      'organism_scientific_name' and 'cluster_representative'.
    """
    # Convert to categorical types for efficiency
    df["organism_scientific_name"] = df["organism_scientific_name"].astype("category")
    df["cluster_representative"] = df["cluster_representative"].astype("category")

    # Shuffle the dataframe to randomize with a fixed seed
    df = df.sample(frac=1, random_state=seed).reset_index(drop=True)

    # Use groupby and head to select the first occurrence after shuffling for each group
    return df.groupby(
        ["organism_scientific_name", "cluster_representative"], observed=True
    ).head(1)


# Defining file paths
labels_path = "data/processed_data/af2/labels.csv"
labels_non_redundant_path = "data/raw_data/af2_structures_non_redundant.csv"
data_exploration_path = "analysis/data_exploration/af2/nonredundant/"


# Reading in label data
labels = pd.read_csv(labels_path)

# Use the optimized function with a seed
labels_non_redundant = select_random_by_organism_and_cluster_optimized(labels, seed=42)

# Saving data set
labels_non_redundant.to_csv(labels_non_redundant_path, index=False)


# Data exploration

# Count the number of rows for each organism before and after
initial_counts = labels["organism_scientific_name"].value_counts()
non_redundant_counts = labels_non_redundant["organism_scientific_name"].value_counts()

# Count unique cluster representatives by organism
unique_clusters_initial = labels.groupby("organism_scientific_name", observed=True)[
    "cluster_representative"
].nunique()
unique_clusters_non_redundant = labels_non_redundant.groupby(
    "organism_scientific_name", observed=True
)["cluster_representative"].nunique()

# Dataframe to hold the summary results
summary_df = pd.DataFrame(
    {
        "Initial_Row_Counts": initial_counts,
        "Non_Redundant_Row_Counts": non_redundant_counts,
        "Initial_Unique_Clusters": unique_clusters_initial,
        "Non_Redundant_Unique_Clusters": unique_clusters_non_redundant,
    }
)

# Outputting summary df
summary_df.to_csv(data_exploration_path + "non_redundant_protein_check.csv")
