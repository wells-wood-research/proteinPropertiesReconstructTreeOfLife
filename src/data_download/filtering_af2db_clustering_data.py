# 0. Importing packages---------------------------------------
import pandas as pd

# 1. Defining variables---------------------------------------

# Defining the path for the raw data
raw_data_path = "data/raw_data/"

# Defining the path for the processed data
processed_data_path = "data/processed_data/af2/"

# Defining the file path for the af2db clusters file
af2db_clusters_file = raw_data_path + "5-allmembers-repId-entryId-cluFlag-taxId.tsv"

# Defining the file path to the processed af2 destress label data
processed_af2_destress_label_data_path = processed_data_path + "labels.csv"

# Defining the output path
output_path = raw_data_path + "filtered_af2db_clusters_data.csv"

# 2. Reading in processed af2 data--------------------------------------

# Reading in data
processed_af2_destress_label_data = pd.read_csv(processed_af2_destress_label_data_path)

# Extracting the uniprot id from the design name in the DE-STRESS data so that we can join on the AF2DB and uniprot data
processed_af2_destress_label_data["uniprot_id"] = (
    processed_af2_destress_label_data["design_name"].str.split("-").str[1]
)

# Extracting a list of the uniprot ids
uniprot_id_list = processed_af2_destress_label_data["uniprot_id"].to_list()

# Initialize a CSV writer if needed, or handle the first chunk differently
first_chunk = True

# Defining chunk count
chunk_count = 0

# Read the large TSV file in chunks
for chunk in pd.read_csv(
    af2db_clusters_file,
    sep="\t",
    chunksize=500000,
    header=None,
):  # Adjust chunksize based on your system's memory
    # Filter the chunk
    filtered_chunk = chunk[chunk.iloc[:, 1].isin(uniprot_id_list)]

    # Process or save the filtered chunk
    # Append to a CSV file without header after first chunk
    filtered_chunk.to_csv(output_path, mode="a", header=first_chunk, index=False)

    chunk_count = chunk_count + 1

    print(chunk_count)
