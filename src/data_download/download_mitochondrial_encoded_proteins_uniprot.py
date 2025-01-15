import requests
import pandas as pd
from dataclasses import dataclass
from typing import Optional, List


def get_response_json(url, params):
    """
    Makes a GET request to the specified URL with the given parameters and returns the JSON response.
    """
    response = requests.get(url, params=params)
    response.raise_for_status()  # This will raise an error for HTTP codes 400/500, which should be handled in calling function
    return response.json()


def get_uniprot_results(query):
    """
    Fetches results from UniProt API for given query parameters.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "json",
        "fields": "accession, protein_name, organism_name, organelle",
    }
    response = requests.get(base_url, params=params)
    response.raise_for_status()

    results = response.json().get("results", [])

    if not results:
        return {}

    return results


# query = 'proteomecomponent:"Mitochondrion MT"'
query = 'organelle: "Mitochondrion"'
results = get_uniprot_results(query)

data = []
for entry in results:
    primary_accession = entry["primaryAccession"]
    scientific_name = entry["organism"]["scientificName"]
    protein_name = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]

    # Identify if the gene encoding type includes 'Mitochondrion'
    gene_encoding_types = [loc["geneEncodingType"] for loc in entry["geneLocations"]]
    # is_mitochondrial = "Mitochondrion" in gene_encoding_types

    # Append the extracted data as a tuple to the list named 'data'
    data.append((primary_accession, protein_name, scientific_name, gene_encoding_types))

# Create a DataFrame
df = pd.DataFrame(
    data,
    columns=[
        "Primary Accession",
        "Protein Name",
        "Organism",
        "Gene Encoded on",
    ],
)

# Defining output path
output_path = "treeoflife/data/raw_data/"

# Saving to a CSV file
df.to_csv(output_path + "uniprot_data_mitochondrial_encoded_proteins.csv", index=False)
