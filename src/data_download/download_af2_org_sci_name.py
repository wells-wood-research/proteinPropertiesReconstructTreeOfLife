# 0. Importing packages and defining custom functions/classes
from dataclasses import dataclass
from typing import Any, Dict, Optional, List
import requests
import pandas as pd
import multiprocessing as mp


@dataclass
class AF2DBResults:
    uniprot_id: Optional[str]
    organism_scientific_name: Optional[str]
    uniprot_description: Optional[str]


def get_response_json(url, params):
    """
    Makes a GET request to the specified URL with the given parameters and returns the JSON response.
    """
    response = requests.get(url, params=params)
    response.raise_for_status()
    return response.json()


def get_af2db_results(uniprot_id):
    """
    Parses results from UniProt query.
    """

    base_url = "https://alphafold.ebi.ac.uk/api/prediction/"
    uniprot_id_url = base_url + uniprot_id

    response = requests.get(uniprot_id_url)
    response.raise_for_status()

    uniprot_description = response.json()[0].get("uniprotDescription", [])
    organism_scientific_name = response.json()[0].get("organismScientificName", [])

    results = AF2DBResults(
        uniprot_id=uniprot_id,
        organism_scientific_name=organism_scientific_name,
        uniprot_description=uniprot_description,
    )

    results_dict = results.__dict__

    return results_dict


def download_af2db_data(
    uniprot_id_list: List[str], output_path: str, num_processes: int
):

    # Create a multiprocessing Pool
    with mp.Pool(processes=num_processes) as pool:
        # Use map to distribute the workload
        results_list = pool.map(get_af2db_results, uniprot_id_list)

    # Convert results dictionary into a DataFrame
    results_df = pd.DataFrame(results_list)
    results_df.to_csv(output_path + "af2db_org_uniprot_desc_results.csv", index=False)

    return results_df


# 1. Defining variables---------------------------------------------------------------------------------

# Defining file path for raw DE-STRESS data
raw_destress_data_af2_path = "data/raw_data/destress_data_af2.csv"

# Defining output path
output_path = "data/raw_data/"

# 2. Loading in data and extracting list of uniprot ids--------------------------------------------------

# Loading in DE-STRESS data
raw_destress_data_af2 = pd.read_csv(raw_destress_data_af2_path)

# Splitting out uniprot id from file name
af2_file_names = raw_destress_data_af2["design_name"].str.split("-").str[1]

# Extracting this as a list
af2_uniprot_id_list = af2_file_names.to_list()

af2_uniprot_id_list = af2_uniprot_id_list[0:100]

# Downloading Uniprot data
results_df = download_af2db_data(
    uniprot_id_list=af2_uniprot_id_list,
    output_path=output_path,
    num_processes=4,
)

print(results_df)
