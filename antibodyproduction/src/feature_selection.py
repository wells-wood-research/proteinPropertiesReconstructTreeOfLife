# This script uses different feature selection methods
# on the Fleishman DE-STRESS metrics.

# 0. Importing packages----------------------------------------------------------
from feature_selection_tools import *

# 1. Defining variables------------------------------------------------------------

# Defining the predictor variable
predictor_var = "expression_bin"

# Defining the feature selection path
feature_selection_path = "antibodyproduction/feature_selection/"

# Defining the scaling methods list
scaling_method_list = ["standard", "robust", "minmax"]

# Defining a composition metrics included flag
comp_flag_list = ["comp", "no_comp"]

# Setting random seed
np.random.seed(42)

for scaling_method in scaling_method_list:
    for comp_flag in comp_flag_list:
        train_data_path = (
            "antibodyproduction/data/processed_data/"
            + scaling_method
            + "/"
            + comp_flag
            + "/"
        )
        feature_selection_scaler_path = (
            feature_selection_path + scaling_method + "/" + comp_flag + "/"
        )

        # 2. Reading in data-------------------------------------------------------------

        # Loading in X_train
        X_train_scaled = pd.read_csv(train_data_path + "X_train_scaled.csv")

        # Loading in y_train
        y_train = pd.read_csv(train_data_path + "y_train.csv")

        # Filtering the y_train to just include the design_name and predictor variable
        y_train = y_train[["design_name", predictor_var]]

        # Extracting the predictor variable on its own
        y_train_pred_var = y_train[predictor_var]

        # 3. Feature selection------------------------------------------------------------

        # Defining the full set of features
        columns_dict = {
            "feature_names": X_train_scaled.columns,
        }
        selected_features_all = pd.DataFrame(columns_dict)
        print(selected_features_all)

        # Using mutual information for feature selection
        selected_features_mi = feature_select_mi(X_train_scaled, y_train_pred_var)

        print(selected_features_mi)

        # Using a random forrest to feature select
        selected_features_rf = feature_select_rf(X_train_scaled, y_train_pred_var)

        print(selected_features_rf)

        # Outputting these list of selected features
        selected_features_all.to_csv(
            feature_selection_scaler_path + "selected_features_all.csv", index=False
        )
        selected_features_mi.to_csv(
            feature_selection_scaler_path + "selected_features_mi.csv", index=False
        )
        selected_features_rf.to_csv(
            feature_selection_scaler_path + "selected_features_rf.csv", index=False
        )
