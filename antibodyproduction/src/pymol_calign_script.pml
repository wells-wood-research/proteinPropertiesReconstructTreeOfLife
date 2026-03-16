# Change to your directory containing the PDB files
cd /home/michael/GitRepos/proteinPropertiesReconstructTreeOfLife/antibodyproduction/data/AF2_structural_models/batch2/

# Start Python block for loading and processing PDB files
python
import os

# Load all PDB files that contain '1ins' in the filename, and include 'rank_001'
pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb') and ('1ins' in f or '5ins' in f) and 'rank_001' in f]

# Load selected PDB files
for pdb in pdb_files:
    cmd.load(pdb)
python end

# Get the list of all loaded objects
all_objects = cmd.get_object_list('all')

# Align all models to the first one
python
for i in range(1, len(all_objects)):
    cmd.align(all_objects[i], all_objects[0])
python end
