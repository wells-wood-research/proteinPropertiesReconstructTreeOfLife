# import requests


# def get_af2db_structures(uniprot_id):
#     """
#     Parses results from UniProt query.
#     """

#     base_url = "https://alphafold.ebi.ac.uk/files/"
#     uniprot_id_url = base_url + "AF-" + uniprot_id + "-F1-model_v4.pdb"

#     response = requests.get(uniprot_id_url)
#     response.raise_for_status()


# get_af2db_structures(uniprot_id="P00403")


import requests
import os
import pandas as pd


def get_af2db_structures(
    uniprot_id,
    # download_dir="analysis/pca_single_proteins/af2/standard/Cytochrome b/af2db_files",
    download_dir="analysis/pca_single_af2db_cluster/af2/standard/A0A7Y5N281/af2db_files",
):
    """
    Downloads PDB files from AlphaFold Database based on UniProt ID.
    """
    base_url = "https://alphafold.ebi.ac.uk/files/"
    pdb_filename = f"{uniprot_id}.pdb"
    uniprot_id_url = base_url + pdb_filename

    try:
        response = requests.get(uniprot_id_url)
        response.raise_for_status()  # Checks if the request returned an HTTP error

        # Create the directory if it does not exist
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)

        # Construct the full path where the file will be saved
        pdb_path = os.path.join(download_dir, pdb_filename)

        # Save the PDB file
        with open(pdb_path, "wb") as f:
            f.write(response.content)

        print(f"Downloaded {pdb_filename} to {download_dir}")

    except requests.HTTPError as e:
        print(
            f"Failed to download structure for {uniprot_id}: HTTP Error {e.response.status_code}"
        )
    except Exception as e:
        print(f"An error occurred for {uniprot_id}: {e}")


def process_structures_from_csv(csv_file):
    """
    Reads UniProt IDs from a CSV file and downloads corresponding structures.
    """
    # Load UniProt IDs from a CSV file into DataFrame
    df = pd.read_csv(csv_file)

    # Assuming the column containing UniProt IDs is named 'UniProt ID'
    for uniprot_id in df["design_name"]:
        get_af2db_structures(uniprot_id)


if __name__ == "__main__":
    # Path to your CSV file containing UniProt IDs
    csv_file_path = "analysis/pca_single_af2db_cluster/af2/standard/A0A7Y5N281/pca_transformed_data.csv"  # Make sure to provide the correct path here
    process_structures_from_csv(csv_file_path)
