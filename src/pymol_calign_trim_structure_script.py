# Start Python block for loading and processing PDB files
import os
import random
from pymol import cmd, stored


def align_and_trim_start(threshold=2.0):
    """
    Align multiple structures to a reference and trim all to well-aligned residues starting from the sequence.

    Parameters:
        structures (list): List of structure identifiers in PyMOL.
        reference (str): Reference structure name.
        start_residue (int): Start residue number for the segment of interest.
        end_residue (int): End residue number for the segment of interest.
        threshold (float): Threshold for considering residues well-aligned.
    """
    alignments = {}

    # Get the list of all loaded objects
    all_objects = cmd.get_object_list("all")

    reference = all_objects[0]

    for structure in range(1, len(all_objects)):
        alignments[structure] = []
        cmd.align(all_objects[structure], all_objects[0], object=f"align_{structure}")

        # # Align each structure to the reference within specified residue range and create an alignment selection
        # for structure in structures:
        #     alignments[structure] = []
        #     cmd.align(
        #         f"{structure} and resi {start_residue}-{end_residue}",
        #         f"{reference} and resi {start_residue}-{end_residue}",
        #         object=f"align_{structure}",
        #     )

        # Get alignment raw data
        aligned_atoms = cmd.get_raw_alignment(f"align_{structure}")
        for pair in aligned_atoms:
            # Extract the model and index from each tuple
            reference_model, target_idx = pair[0]
            structure_model, mobile_idx = pair[1]

            dist = cmd.get_distance(
                f"{reference_model} and index {target_idx}",
                f"{structure_model} and index {mobile_idx}",
            )

            if dist <= threshold:
                target_resi = (
                    cmd.get_model(f"{reference} and index {target_idx}").atom[0].resi
                )
                mobile_resi = (
                    cmd.get_model(f"{all_objects[structure]} and index {mobile_idx}")
                    .atom[0]
                    .resi
                )
                # # Check if the residues are within the range
                # if (
                #     int(target_resi) >= start_residue
                #     and int(target_resi) <= end_residue
                # ):
                alignments[structure].append((target_resi, mobile_resi))

    # Intersection of well-aligned residues within range
    common_residues = set(alignments[all_objects[0]])
    for structure in all_objects:
        common_residues.intersection_update(set(alignments[structure]))

    if not common_residues:
        print("No common well-aligned residues found within the specified range.")
        return

    # Create the trimmed structures and output their sequences
    intersection_selection = " or ".join(
        [
            f"{structure} and resi {pair[1]}"
            for structure in all_objects
            for pair in common_residues
            # if int(pair[1]) >= start_residue and int(pair[1]) <= end_residue
        ]
    )
    cmd.select("well_aligned_range", intersection_selection)
    for structure in all_objects:
        cmd.create(f"trimmed_{structure}", f"{structure} and well_aligned_range")
        print(f"Sequence for {structure}:")
        print(cmd.get_fastastr(f"trimmed_{structure}"))

    # Process the reference sequence last
    cmd.create(f"trimmed_{reference}", f"{reference} and well_aligned_range")
    print(f"Sequence for {reference}:")
    print(cmd.get_fastastr(f"trimmed_{reference}"))


# # Example usage
# structures = [
#     "AF-A0A1D8PQH5-F1-model_v4",
#     "AF-Q4D9I5-F1-model_v4",
#     "AF-P41979-F1-model_v4",
#     "AF-Q57Z58-F1-model_v4",
#     "AF-Q9LYK8-F1-model_v4",
#     "AF-P41980-F1-model_v4",
#     "AF-P09671-F1-model_v4",
#     "AF-J9EJI4-F1-model_v4",
#     "AF-C0NJR1-F1-model_v4",
#     "AF-A0A175W4W0-F1-model_v4",
#     "AF-P31161-F1-model_v4",
#     "AF-P43019-F1-model_v4",
#     "AF-P28759-F1-model_v4",
#     "AF-A0A175W7X9-F1-model_v4",
#     "AF-A0A077Z4X9-F1-model_v4",
#     "AF-A0A0D2GY51-F1-model_v4",
#     "AF-B4F925-F1-model_v4",
#     "AF-P0A2F4-F1-model_v4",
#     "AF-A0A0N4UAL0-F1-model_v4",
#     "AF-P07895-F1-model_v4",
#     "AF-A0A0H3GLE8-F1-model_v4",
#     "AF-P00448-F1-model_v4",
#     "AF-Q00637-F1-model_v4",
#     "AF-P41978-F1-model_v4",
#     "AF-A0A0D2G8Z9-F1-model_v4",
#     "AF-I1LCI3-F1-model_v4",
#     "AF-P41977-F1-model_v4",
#     "AF-U7Q995-F1-model_v4",
#     "AF-P0A0J3-F1-model_v4",
#     "AF-A0A1D6FNG1-F1-model_v4",
#     "AF-A4HTI0-F1-model_v4",
#     "AF-P53652-F1-model_v4",
#     "AF-Q9UQX0-F1-model_v4",
#     "AF-P13367-F1-model_v4",
#     "AF-P43725-F1-model_v4",
#     "AF-Q6P980-F1-model_v4",
#     "AF-Q32A73-F1-model_v4",
#     "AF-P04179-F1-model_v4",
#     "AF-P00447-F1-model_v4",
#     "AF-C0NJ49-F1-model_v4",
#     "AF-C1GW30-F1-model_v4",
#     "AF-A0A1C1D0I2-F1-model_v4",
#     "AF-P09233-F1-model_v4",
#     "AF-K0EMC4-F1-model_v4",
#     "AF-P43312-F1-model_v4",
#     "AF-Q55BJ9-F1-model_v4",
#     "AF-G4VSQ5-F1-model_v4",
#     "AF-O74379-F1-model_v4",
#     "AF-Q4DCQ3-F1-model_v4",
#     "AF-P9WGE7-F1-model_v4",
#     "AF-A0A0K0EMH0-F1-model_v4",
#     "AF-Q5A8Z4-F1-model_v4",
#     "AF-A0A1C1D000-F1-model_v4",
#     "AF-Q2G261-F1-model_v4",
#     "AF-P41981-F1-model_v4",
# ]  # list your structure identifiers here
# reference = "AF-Q0PBW9-F1-model_v4"
# start_residue, end_residue = (
#     1,
#     50,
# )
# Define the range of residues at the start of the sequence
# Change to your directory containing the PDB files
os.chdir(
    "/home/mstam/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/pca_single_af2db_cluster/af2/standard/A0A0G9LKG2/af2_structures/"
)

# Load all PDB files
pdb_files = [f for f in os.listdir(".") if f.endswith(".pdb")]

# Select 20 random PDB files from the list
selected_pdb_files = random.sample(pdb_files, 3) if len(pdb_files) >= 3 else pdb_files

# Load selected PDB files into your molecular viewer
for pdb_file in pdb_files:
    cmd.load(pdb_file)

align_and_trim_start()
