# Change to your directory containing the PDB files
cd /home/mstam/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/pca_single_af2db_cluster/af2/standard/A0A0G9LKG2/af2_structures/

# Start Python block for loading and processing PDB files
python
import os
import random

# Load all PDB files 
pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]

# Select 20 random PDB files from the list
# selected_pdb_files = random.sample(pdb_files, 20) if len(pdb_files) >= 20 else pdb_files

# Load selected PDB files into your molecular viewer
for pdb_file in pdb_files:
    cmd.load(pdb_file)
python end

# Get the list of all loaded objects
all_objects = cmd.get_object_list('all')

# Align all models to the first one
python
for i in range(1, len(all_objects)):
    cmd.align(all_objects[i], all_objects[0])
python end
