import pandas as pd
import os
import glob


def pdb2df(pdbFile):
    columns = [
        "ATOM",
        "ATOM_ID",
        "ATOM_NAME",
        "RES_NAME",
        "CHAIN_ID",
        "RES_SEQ",
        "X",
        "Y",
        "Z",
        "OCCUPANCY",
        "TEMP_FACTOR",
        "ELEMENT",
    ]
    data = []
    with open(pdbFile, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_type = line[0:6].strip()
                atom_id = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                if chain_id == "":
                    chain_id = None
                res_seq = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occupancy = float(line[54:60].strip())
                temp_factor = float(line[60:66].strip())
                element = line[76:78].strip()

                data.append(
                    [
                        atom_type,
                        atom_id,
                        atom_name,
                        res_name,
                        chain_id,
                        res_seq,
                        x,
                        y,
                        z,
                        occupancy,
                        temp_factor,
                        element,
                    ]
                )

    return pd.DataFrame(data, columns=columns)


def get_mean_af2_pLDDT(pdb_file_path):
    pdbdf = pdb2df(pdb_file_path)
    mean_residue_pLDDT = pdbdf[pdbdf["ATOM_NAME"] == "CA"]["TEMP_FACTOR"].mean()
    return mean_residue_pLDDT


def process_folder(input_folder, output_csv):
    columns = ["Filename", "Mean_PLDDT"]
    result_data = []

    count = 1
    os.chdir(input_folder)

    full_list = len(glob.glob("*.pdb"))

    for filename in glob.glob("*.pdb"):
        print("File no. " + str(count) + "/ " + str(full_list))
        pdb_file_path = os.path.join(input_folder, filename)
        mean_plddt = get_mean_af2_pLDDT(pdb_file_path)
        result_data.append([filename, mean_plddt])
        count = count + 1

    result_df = pd.DataFrame(result_data, columns=columns)
    result_df.to_csv(output_csv, index=False)


# Replace 'input_folder' and 'output.csv' with your actual folder path and output file name
process_folder(
    "/home/michael/GitRepos/illuminating-protein-structural-universe/pdb_af2_structural_analysis/data/raw_data/af2/designed_pdb_results/",
    "plddt_scores_designed_pdbs.csv",
)
