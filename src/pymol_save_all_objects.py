# Python script to run within PyMOL
from pymol import cmd, stored

# Retrieve all objects in the session
objects = cmd.get_names("all")

for obj in objects:
    filename = f"/home/mstam/GitRepos/proteinPropertiesReconstructTreeOfLife/data/raw_data/superoxide_dismutase_trunc/{obj}.pdb"  # Create a filename from each object name
    cmd.save(
        filename, obj
    )  # Save each object as a PDB file with its object name as filename

print("All objects have been exported.")
