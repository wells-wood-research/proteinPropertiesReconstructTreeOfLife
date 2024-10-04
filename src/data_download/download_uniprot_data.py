# 0. Importing packages and defining custom functions/classes
from dataclasses import dataclass, field
from typing import Any, Dict, Generic, Optional, Tuple, TypeVar, List
import requests
import pandas as pd
import multiprocessing as mp


@dataclass
class UniprotResults:
    uniprot_id: Optional[str]
    host_organism: Optional[str]
    organism_class: Optional[str]
    go_codes: Optional[str]
    subcellular_location: Optional[str]


def get_response_json(url, params):
    """
    Makes a GET request to the specified URL with the given parameters and returns the JSON response.
    """
    response = requests.get(url, params=params)
    response.raise_for_status()
    return response.json()


def extract_subcellular_locations(comments):
    """
    Extracts subcellular locations from comments.
    """
    subcellular_locations = []
    for comment in comments:
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            locations = comment.get("subcellularLocations", [])
            subcellular_locations += [loc["location"]["value"] for loc in locations]
    return ":".join(subcellular_locations) if subcellular_locations else ""


def extract_go_codes(cross_references):
    """
    Extracts GO codes from cross references.
    """
    return ":".join(
        ref["id"] for ref in cross_references if ref.get("database") == "GO"
    )


def get_uniprot_results(uniprot_id):
    """
    Parses results from UniProt query.
    """

    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": uniprot_id,
        "format": "json",
    }

    response = requests.get(base_url, params=params)
    response.raise_for_status()

    results = response.json().get("results", [])

    if not results:
        return {}

    result = results[0]
    uniprot_id = result.get("primaryAccession", "")
    organism_info = result.get("organism", {})
    host_organism = organism_info.get("scientificName", "")
    organism_class = (
        organism_info.get("lineage", [])[0] if organism_info.get("lineage", []) else ""
    )

    cross_references = result.get("uniProtKBCrossReferences", [])
    go_codes_str = extract_go_codes(cross_references)

    comments = result.get("comments", [])
    subcellular_locations_str = extract_subcellular_locations(comments)

    results = UniprotResults(
        uniprot_id,
        host_organism,
        organism_class,
        go_codes_str,
        subcellular_locations_str,
    )

    results_dict = results.__dict__

    return results_dict


# def get_uniprot_entry(uniprot_id: str) -> UniprotResults:
#     """
#     Retrieve the UniProt host organism (scientific name) and subcellular location for a UniProt ID

#     Parameters
#     ----------
#     uniprot_id : str
#         The uniprot ID to query.

#     Returns
#     -------
#     UniprotResults
#         A dataclass returning the uniprot id, host organism name and subcellular location
#     """
#     base_url = "https://rest.uniprot.org/uniprotkb/search"
#     params = {
#         "query": uniprot_id,
#         "format": "json",
#     }

#     response = requests.get(base_url, params=params)
#     response.raise_for_status()

#     results = response.json().get("results", [])
#     if results:
#         if "primaryAccession" in results[0]:
#             uniprot_id = results[0]["primaryAccession"]
#         else:
#             uniprot_id = ""
#         if "organism" in results[0]:
#             host_organism = results[0]["organism"]["scientificName"]
#             organism_class = results[0]["organism"]["lineage"][0]

#         else:
#             host_organism = ""
#             organism_class = ""

#         cross_references = results[0].get("uniProtKBCrossReferences", [])

#         go_codes = []
#         for ref in cross_references:
#             if ref["database"] == "GO":
#                 go_codes.append(ref["id"])

#         go_codes_str = ":".join(go_codes)

#         subcellular_locations = []

#         if "comments" in results[0]:
#             comments = results[0]["comments"]
#             for comment in comments:
#                 if comment["commentType"] == "SUBCELLULAR LOCATION":
#                     if "subcellularLocations" in comment:
#                         for subcellular_location in comment["subcellularLocations"]:
#                             subcellular_locations.append(
#                                 subcellular_location["location"]["value"]
#                             )

#             subceullar_locations_str = ":".join(subcellular_locations)

#         else:
#             subceullar_locations_str = ""

#         results = UniprotResults(
#             uniprot_id,
#             host_organism,
#             organism_class,
#             go_codes_str,
#             subceullar_locations_str,
#         )

#         results_dict = results.__dict__

#         return results_dict


# def get_uniprot_dict(uniprot_id: str) -> dict:

#     try:
#         results = get_uniprot_entry(uniprot_id)
#         print(results)
#         return {uniprot_id: results.__dict__}
#     except Exception:
#         return {uniprot_id: {}}


def download_uniprot_data(
    uniprot_id_list: List[str], output_path: str, num_processes: int
):

    # Create a multiprocessing Pool
    with mp.Pool(processes=num_processes) as pool:
        # Use map to distribute the workload
        results_list = pool.map(get_uniprot_results, uniprot_id_list)

    # Convert results dictionary into a DataFrame
    results_df = pd.DataFrame(results_list)
    results_df.to_csv(output_path + "uniprot_results_org_subcellloc.csv", index=False)

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
results_df = download_uniprot_data(
    uniprot_id_list=af2_uniprot_id_list,
    output_path=output_path,
    num_processes=4,
)

print(results_df)
