# 0. Importing packages and helper functions---------------------------------------------
import numpy as np
import pandas as pd

# 1. Defining variables------------------------------------------------------------------

# Defining path to raw uniprot data
raw_uniprot_data_path = "data/raw_data/uniprot_results_org_subcellloc_uniprotkb.csv"

# Defining path to raw af2db data
raw_af2db_data_path = "data/raw_data/af2db_org_uniprot_desc_results.csv"

# Defining the output paths
data_output_path = "data/processed_data/"
analysis_output_path = "analysis/data_exploration/uniprot/"

# 2. Reading in data----------------------------------------------------------------------

# Reading in raw uniprot data
raw_uniprot_data = pd.read_csv(raw_uniprot_data_path)

# Reading in raw af2db data
raw_af2db_data = pd.read_csv(raw_af2db_data_path)

print(raw_uniprot_data)
print(raw_af2db_data)

# 3. Processing data----------------------------------------------------------------------

# First removing duplicates from uniprot data and af2db data
raw_uniprot_data = raw_uniprot_data.drop_duplicates().reset_index(drop=True)
raw_af2db_data = raw_af2db_data.drop_duplicates().reset_index(drop=True)

print(raw_uniprot_data.shape)
print(raw_af2db_data.shape)

# Joining these data sets together by uniprot id
af2db_uniprot_data = pd.merge(
    raw_af2db_data,
    raw_uniprot_data[["uniprot_id", "host_organism", "subcellular_location"]],
    on="uniprot_id",
    how="left",
)

print(af2db_uniprot_data.shape)

# Setting the missing org sci names from af2db to the values from host_organism from uniprot. These are all homo sapiens for some reason.
af2db_uniprot_data.loc[
    af2db_uniprot_data["organism_scientific_name"].isna(), "organism_scientific_name"
] = af2db_uniprot_data["host_organism"]

print(af2db_uniprot_data)


# Test for mismatching org names
test = af2db_uniprot_data[
    af2db_uniprot_data["organism_scientific_name"]
    != af2db_uniprot_data["host_organism"]
].reset_index(drop=True)

test.to_csv("mismatch_organism_test.csv")

# Etracting the org name and removing anything else including in this field (e.g strain etc.)
af2db_uniprot_data["organism_scientific_name"] = (
    af2db_uniprot_data["organism_scientific_name"]
    .str.split(" ")
    .str[0:2]
    .apply(lambda x: " ".join(x) if x is not np.nan else "Unknown")
)


# Data exploration
af2db_uniprot_data.value_counts("uniprot_description").to_csv(
    analysis_output_path + "af2db_uniprot_data_uniprot_desc_count.csv"
)
af2db_uniprot_data.value_counts("subcellular_location").to_csv(
    analysis_output_path + "af2db_uniprot_data_subcell_loc_count.csv"
)
af2db_uniprot_data.value_counts("organism_scientific_name").to_csv(
    analysis_output_path + "af2db_uniprot_data_org_name_count.csv"
)

# Identifying the cytoplasm subceullar location
af2db_uniprot_data["comments_cytoplasm_flag"] = af2db_uniprot_data[
    "subcellular_location"
].str.contains("Cytoplasm")

# Filtering to cytoplasm only
af2db_uniprot_data_cytoplasm = af2db_uniprot_data[
    af2db_uniprot_data["comments_cytoplasm_flag"] == 1
].reset_index(drop=True)

# Distribution of org name for cytoplasm proteins
af2db_uniprot_data_cytoplasm.value_counts("organism_scientific_name").to_csv(
    analysis_output_path + "af2db_uniprot_data_org_name_cyto_count.csv"
)

# Filtering to Protein kinase domain-containing protein only
af2db_uniprot_data_kinase = af2db_uniprot_data[
    af2db_uniprot_data["uniprot_description"]
    == "Protein kinase domain-containing protein"
].reset_index(drop=True)

# Distribution of org name for kinase proteins
af2db_uniprot_data_kinase.value_counts("organism_scientific_name").to_csv(
    analysis_output_path + "af2db_uniprot_data_org_name_kinase_count.csv"
)

# Identifying protein deescriptions that are present across all organisms

# Group by protein name and aggregate unique organisms
protein_counts = af2db_uniprot_data.groupby("uniprot_description")[
    "organism_scientific_name"
].nunique()

# Filter for proteins that appear in greater than 40 organisms
proteins_in_all_organisms = protein_counts[protein_counts >= 40]
proteins_in_all_organisms.sort_values(inplace=True, ascending=False)
proteins_in_all_organisms.to_csv(
    analysis_output_path + "uniprot_desc_org_count_gt_40.csv"
)


# Filtering to tRNA (guanine-N(7)-)-methyltransferase only
af2db_uniprot_data_trna = af2db_uniprot_data[
    af2db_uniprot_data["uniprot_description"]
    == "tRNA (guanine-N(7)-)-methyltransferase"
].reset_index(drop=True)

# Distribution of org name for kinase proteins
af2db_uniprot_data_trna.value_counts("organism_scientific_name").to_csv(
    analysis_output_path + "af2db_uniprot_data_org_name_trna_count.csv"
)


# Removing host organism field
af2db_uniprot_data.drop(["host_organism"], inplace=True, axis=1)


# Outputting data
af2db_uniprot_data.to_csv(
    data_output_path + "processed_af2db_uniprot_data.csv", index=False
)


# af2db_uniprot_data = raw_uniprot_data[
#     ~raw_uniprot_data["organism_scientific_name"].isin(
#         ["Unknown", "Metarhizium anisopliae", "Neovison vison"]
#     )
# ].reset_index(drop=True)

# print(raw_uniprot_data)

# # Define the specific GO code
# cytoplasm_go_code = "GO:0005737"


# raw_uniprot_data["go_codes"] = raw_uniprot_data["go_codes"].str.replace(":GO", ";GO")


# raw_uniprot_data["go_cytoplasm_flag"] = raw_uniprot_data["go_codes"].apply(
#     lambda x: 1 if cytoplasm_go_code in (x.split(";") if pd.notna(x) else []) else 0
# )


# # raw_uniprot_data["comments_cytoplasm_flag"]

# cytoplasm_proteins = raw_uniprot_data[
#     (raw_uniprot_data["go_cytoplasm_flag"] == 1)
#     | (raw_uniprot_data["comments_cytoplasm_flag"] == 1)
# ].reset_index(drop=True)

# print(cytoplasm_proteins)

# print(cytoplasm_proteins.value_counts("organism_scientific_name"))
