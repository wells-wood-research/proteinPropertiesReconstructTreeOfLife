cd /home/mstam/GitRepos/proteinPropertiesReconstructTreeOfLife/analysis/pca_single_af2db_cluster/af2/standard/A0A0G9LKG2/af2_structures/

python
import os

pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]

for pdb_file in pdb_files:
    cmd.load(pdb_file)
python end

python
all_objects = cmd.get_object_list('all')

for i in range(1, len(all_objects)):
    cmd.align(all_objects[i], all_objects[0])
python end

hide everything, resi 1-70 and b < 70
