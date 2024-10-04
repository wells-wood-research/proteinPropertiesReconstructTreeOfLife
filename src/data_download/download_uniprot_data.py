# 0. Importing packages and defining custom functions/classes
from dataclasses import dataclass, field
from typing import Any, Dict, Generic, Optional, Tuple, TypeVar, List
import requests
import pandas as pd


@dataclass
class UniprotResults:
    uniprot_id: str
    host_organism: Optional[str]
    organism_class: Optional[str]
    go_codes: Optional[str]
    subcellular_location: Optional[str]


def get_uniprot_entry(uniprot_id: str) -> UniprotResults:
    """
    Retrieve the UniProt host organism (scientific name) and subcellular location for a UniProt ID

    Parameters
    ----------
    uniprot_id : str
        The uniprot ID to query.

    Returns
    -------
    UniprotResults
        A dataclass returning the uniprot id, host organism name and subcellular location
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": uniprot_id,
        "format": "json",
    }

    response = requests.get(base_url, params=params)
    response.raise_for_status()

    results = response.json().get("results", [])
    if results:
        uniprot_id = results[0]["primaryAccession"]
        if "organism" in results[0]:
            host_organism = results[0]["organism"]["scientificName"]
            organism_class = results[0]["organism"]["lineage"][0]

        else:
            host_organism = ""
            organism_class = ""

        cross_references = results[0].get("uniProtKBCrossReferences", [])

        go_codes = []
        for ref in cross_references:
            if ref["database"] == "GO":
                go_codes.append(ref["id"])

        go_codes_str = ":".join(go_codes)

        subcellular_locations = []

        if "comments" in results[0]:
            comments = results[0]["comments"]
            for comment in comments:
                if comment["commentType"] == "SUBCELLULAR LOCATION":
                    if "subcellularLocations" in comment:
                        for subcellular_location in comment["subcellularLocations"]:
                            subcellular_locations.append(
                                subcellular_location["location"]["value"]
                            )

            subceullar_locations_str = ":".join(subcellular_locations)

        else:
            subceullar_locations_str = ""

        return UniprotResults(
            uniprot_id,
            host_organism,
            organism_class,
            go_codes_str,
            subceullar_locations_str,
        )
    else:
        raise Exception(f"No entry found for uniprot ID: {uniprot_id}")


def download_uniprot_data(uniprot_id_list: List[str], output_path: str):

    results_dict = {}

    for uniprot_id in uniprot_id_list:

        results = get_uniprot_entry(uniprot_id=uniprot_id)
        results_dict[uniprot_id] = results.__dict__

    results_df = pd.DataFrame(results_dict).transpose().reset_index(drop=True)
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

# Downloading Uniprot data
results_df = download_uniprot_data(
    uniprot_id_list=af2_uniprot_id_list, output_path=output_path
)

print(results_df)
