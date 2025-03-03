import csv

# Assuming the data is stored in a CSV file proteins.csv
with open(
    "/home/mstam/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/pca_single_af2db_cluster/af2/standard/A0A0G9LKG2/pca_transformed_data.csv",
    newline="",
) as csvfile, open(
    "/home/mstam/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/pca_single_af2db_cluster/af2/standard/A0A0G9LKG2/A0A0G9LKG2_af2db_cluster_seqs.fasta",
    "w",
) as fastafile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        header = f">{row['design_name']}\n"
        sequence = f"{row['full_sequence']}\n"
        fastafile.write(header)
        fastafile.write(sequence)
