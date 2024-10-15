# This script is the main script which prepares the AF2 and PDB DE-STRESS data
# so that it is ready for downstream analysis.

# 0. Importing packages and helper functions---------------------------------------------
from data_prep_tools import *
import numpy as np


# 1. Defining variables-------------------------------------------------------------------

# Defining the data set list
dataset_list = ["af2", "pdb"]

# Defining the scaling methods list
scaling_method_list = ["standard", "robust", "minmax"]

# Defining the path for the raw data
raw_data_path = "data/raw_data/"

# Defining the path for the processed data
processed_data_path = "data/processed_data/"

# Defining the file path for the AF2 raw destress data
raw_destress_data_af2_path = raw_data_path + "destress_data_af2.csv"

# Defining the file path for the PDB raw destress data
raw_destress_data_pdb_path = raw_data_path + "destress_data_pdb_082024.csv"

# Defining the file path for the af2db uniprot data
af2db_uniprot_data_path = processed_data_path + "processed_af2db_uniprot_data.csv"

# Defining the file path for the plddt scores
af2_plddt_scores_path = raw_data_path + "af2_plddt_scores.csv"

# Defining a list of DE-STRESS metrics which are energy field metrics
energy_field_list = [
    "hydrophobic_fitness",
    "budeff_total",
    "budeff_steric",
    "budeff_desolvation",
    "budeff_charge",
    "evoef2_total",
    "evoef2_ref_total",
    "evoef2_intraR_total",
    "evoef2_interS_total",
    "evoef2_interD_total",
    "rosetta_total",
    "rosetta_fa_atr",
    "rosetta_fa_rep",
    "rosetta_fa_intra_rep",
    "rosetta_fa_elec",
    "rosetta_fa_sol",
    "rosetta_lk_ball_wtd",
    "rosetta_fa_intra_sol_xover4",
    "rosetta_hbond_lr_bb",
    "rosetta_hbond_sr_bb",
    "rosetta_hbond_bb_sc",
    "rosetta_hbond_sc",
    "rosetta_dslf_fa13",
    "rosetta_rama_prepro",
    "rosetta_p_aa_pp",
    "rosetta_fa_dun",
    "rosetta_omega",
    "rosetta_pro_close",
    "rosetta_yhh_planarity",
]

# Defining cols to drop
drop_cols = [
    "ss_prop_alpha_helix",
    "ss_prop_beta_bridge",
    "ss_prop_beta_strand",
    "ss_prop_3_10_helix",
    "ss_prop_pi_helix",
    "ss_prop_hbonded_turn",
    "ss_prop_bend",
    "ss_prop_loop",
    "charge",
    "mass",
    "num_residues",
    "uniprot_id",
    "aggrescan3d_total_value",
    "rosetta_pro_close",
    "rosetta_omega",
    "rosetta_total",
    "rosetta_fa_rep",
    "evoef2_total",
    "evoef2_interS_total",
    "rosetta_rama_prepro",
    "rosetta_p_aa_pp",
    "evoef2_ref_total",
    "Mean_PLDDT",
    "subcellular_location",
]

# Defining the organism groups
organism_animal_list = [
    "Caenorhabditis elegans",
    "Danio rerio",
    "Drosophila melanogaster",
    "Mus musculus",
    "Rattus norvegicus",
    "Homo sapiens",
    "Brugia malayi",
    "Dracunculus medinensis",
    "Onchocerca volvulus",
    "Schistosoma mansoni",
    "Strongyloides stercoralis",
    "Trichuris trichiura",
    "Wuchereria bancrofti",
]

organism_fungi_list = [
    "Candida albicans",
    "Saccharomyces cerevisiae",
    "Cladophialophora carrionii",
    "Fonsecaea pedrosoi",
    "Madurella mycetomatis",
    "Sporothrix schenckii",
    "Ajellomyces capsulatus",
    "Schizosaccharomyces pombe",
    "Paracoccidioides lutzii",
]
organism_bacteria_list = [
    "Escherichia coli",
    "Helicobacter pylori",
    "Campylobacter jejuni",
    "Enterococcus faecium",
    "Klebsiella pneumoniae",
    "Mycobacterium leprae",
    "Mycobacterium tuberculosis",
    "Mycobacterium ulcerans",
    "Neisseria gonorrhoeae",
    "Nocardia brasiliensis",
    "Pseudomonas aeruginosa",
    "Salmonella typhimurium",
    "Shigella dysenteriae",
    "Staphylococcus aureus",
    "Streptococcus pneumoniae",
    "Haemophilus influenzae",
]
organism_plant_list = [
    "Arabidopsis thaliana",
    "Glycine max",
    "Oryza sativa",
    "Zea mays",
]
organism_protozoan_list = [
    "Plasmodium falciparum",
    "Dictyostelium discoideum",
    "Leishmania infantum",
    "Trypanosoma brucei",
    "Trypanosoma cruzi",
]

organism_archaea_list = ["Methanocaldococcus jannaschii"]

organism_other_list = [
    "Other",
    "Unknown",
]

# Defining a dictionary with these organism group lists
organism_group_dict = {
    "animal": organism_animal_list,
    "archaea": organism_archaea_list,
    "bacteria": organism_bacteria_list,
    "fungi": organism_fungi_list,
    "plant": organism_plant_list,
    "protozoan": organism_protozoan_list,
    "other": organism_other_list,
}

# Defining the labels that we are interested in
labels = [
    "design_name",
    "full_sequence",
    "dssp_bin",
    "charge",
    "isoelectric_point",
    "isoelectric_point_bin",
    "rosetta_total",
    "packing_density",
    "packing_density_bin",
    "hydrophobic_fitness",
    "aggrescan3d_avg_value",
    "aggrescan3d_avg_bin",
    "organism_scientific_name",
    "organism_group",
    "organism_group2",
    "uniprot_description",
    "subcellular_location",
    "Mean_PLDDT",
]

# Defining a threshold for the spearman correlation coeffient
# in order to remove highly correlated variables
corr_coeff_threshold = 0.6

# Defining a threshold to remove features that have pretty much the same value
constant_features_threshold = 0.25

# Defining a threshold for removing missing values for af2 data
missing_val_threshold_af2 = 0.05

# Defining a path for the data exploration for af2
data_exploration_af2_path = "analysis/data_exploration/af2/"

# Defining a path for the data output path for af2
processed_data_af2_path = "data/processed_data/af2/"

# Setting a flat to remove low quality af2 models
remove_low_quality_af2_models = True


# 2. Reading in data sets-------------------------------------------------------------------------------

# Reading in raw AF2 DE-STRESS data
raw_destress_data_af2 = pd.read_csv(raw_destress_data_af2_path)

# Reading in af2db uniprot data
af2db_uniprot_data = pd.read_csv(af2db_uniprot_data_path)

# Reading in the plddt scores for the af2 structural models
af2_plddt_scores = pd.read_csv(af2_plddt_scores_path)
af2_plddt_scores["design_name"] = (
    af2_plddt_scores["Filename"].str.replace(".pdb", "").astype(str)
)

# Joining this score onto the af2 structural model data set
raw_destress_data_af2 = raw_destress_data_af2.merge(
    af2_plddt_scores[["design_name", "Mean_PLDDT"]], on="design_name", how="left"
)

if remove_low_quality_af2_models:
    raw_destress_data_af2 = raw_destress_data_af2[
        raw_destress_data_af2["Mean_PLDDT"] >= 70
    ].reset_index(drop=True)
else:
    raw_destress_data_af2 = raw_destress_data_af2


# Reading in raw PDB DE-STRESS data
raw_destress_data_pdb = pd.read_csv(raw_destress_data_pdb_path)

# 3. Processing data sets-------------------------------------------------------------------------------


# AF2

process_af2_data(
    raw_destress_data=raw_destress_data_af2,
    af2db_uniprot_data=af2db_uniprot_data,
    data_exploration_path=data_exploration_af2_path,
    data_output_path=processed_data_af2_path,
    missing_val_threshold=missing_val_threshold_af2,
    organism_group_dict=organism_group_dict,
    energy_field_list=energy_field_list,
    labels=labels,
    drop_cols_list=drop_cols,
    constant_features_threshold=constant_features_threshold,
    scaling_method_list=scaling_method_list,
    corr_coeff_threshold=corr_coeff_threshold,
    remove_low_quality_af2_models=remove_low_quality_af2_models,
)
