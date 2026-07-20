import subprocess
import pandas as pd


def sync_files(server, username, remote_folder, local_folder, file_list):
    # Create a pattern for rsync to include only specified files
    include_patterns = [f"--include={filename}" for filename in file_list]
    include_patterns.append("--exclude=*")  # Exclude all other files

    # Formulate the rsync command
    remote_path = f"{username}@{server}:{remote_folder}/"

    # Build the rsync command
    rsync_command = ["rsync", "-avz", *include_patterns, remote_path, local_folder]

    # Execute the command
    print(f"Syncing selected files from {remote_path} to {local_folder}")

    try:
        subprocess.run(rsync_command, check=True)

    except subprocess.CalledProcessError as e:
        print("ERROR:", e.stderr)


def main():
    server = "glutamate.bio.ed.ac.uk"
    username = "mjstam"
    remote_folder = "/mnt/scratch/alphafold_model_organisms/"
    # local_folder = "/home/mstam/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/pca_single_af2db_cluster/af2/standard/A0A0G9LKG2/af2_structures/"
    local_folder = "/home/mstam/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/pca_single_af2db_cluster/af2/standard/A0A0B8NHG6/af2_structures/"

    pca_data_single_cluster = pd.read_csv(
        "/home/mstam/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/pca_single_af2db_cluster/af2/standard/A0A0B8NHG6/pca_transformed_data.csv",
        # "/home/mstam/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/pca_single_af2db_cluster/af2/standard/A0A0G9LKG2/pca_transformed_data.csv",
    )
    file_names = [
        file_name + ".pdb"
        for file_name in pca_data_single_cluster["design_name"].to_list()
    ]

    # # Add .pdb suffix
    # pdb_filenames = ["AF-" + filename + "-F1-model_v4.pdb" for filename in file_names]

    sync_files(server, username, remote_folder, local_folder, file_names)


if __name__ == "__main__":
    main()
