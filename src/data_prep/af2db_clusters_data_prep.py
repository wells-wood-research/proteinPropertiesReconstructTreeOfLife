# 0. Importing packages---------------------------------------
import pandas as pd

# 1. Defining variables---------------------------------------

# Defining file path for the filtered af2db clusters data
af2db_filtered_clusters_data_path = "data/raw_data/filtered_af2db_clusters_data.csv"

# Defining a path for the processed data
processed_data_path = "data/processed_data/af2/"

# 2. Reading in data--------------------------------------------

# Reading in af2db filtered clusters data
af2db_filtered_clusters_data = pd.read_csv(af2db_filtered_clusters_data_path)
af2db_filtered_clusters_data.columns = [
    "cluster_representative",
    "unuiprot_id",
    "cluster_flag",
    "taxonomic_id",
]

print(af2db_filtered_clusters_data)

# 3. Processing data----------------------------------------------

# Removing incorrect rows
af2db_filtered_clusters_data = af2db_filtered_clusters_data[
    af2db_filtered_clusters_data["cluster_representative"] != "0"
].reset_index(drop=True)

# Removing duplicates
af2db_filtered_clusters_data = (
    af2db_filtered_clusters_data.drop_duplicates().reset_index(drop=True)
)

# Outputting csv file
af2db_filtered_clusters_data.to_csv(
    processed_data_path + "processed_af2db_clusters_data.csv", index=False
)
