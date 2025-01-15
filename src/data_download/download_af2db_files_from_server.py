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
    server = ""
    username = ""
    remote_folder = ""
    local_folder = ""

    pca_data_single_cluster = pd.read_csv(
        "",
        delimiter="\t",
    )
    file_names = pca_data_single_cluster["Entry"].to_list()

    # Add .pdb suffix
    pdb_filenames = ["AF-" + filename + "-F1-model_v4.pdb" for filename in file_names]

    sync_files(server, username, remote_folder, local_folder, pdb_filenames)


if __name__ == "__main__":
    main()
