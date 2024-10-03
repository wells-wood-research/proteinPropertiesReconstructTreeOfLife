# This script is the main script which prepares the AF2 and PDB DE-STRESS data
# so that it is ready for downstream analysis.

# 0. Importing packages and helper functions---------------------------------------------
from data_prep_tools import *
import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler, RobustScaler
import pickle

# 1. Defining variables-------------------------------------------------------------------

# Defining the data set list
dataset_list = ["af2", "pdb"]

# Defining the scaling methods list
scaling_method_list = ["standard", "robust", "minmax"]

# Defining the path for the raw data
raw_data_path = "data/raw_data"

# Defining the file path for the AF2 raw destress data
raw_destress_data_af2_path = raw_data_path + "destress_data_af2.csv"

# Defining the file path for the PDB raw destress data
raw_destress_data_pdb_path = raw_data_path + "destress_data_pdb_082024.csv"

# Defining the file path for the plddt scores
af2_plddt_scores_path = raw_data_path + "af2_plddt_scores.csv"
