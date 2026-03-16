# This script processes the fleishman expression data and structural feature data
# in order to create training and test data sets for machine learning.

# 0. Importing packages------------------------------------------------------------------------
import numpy as np
import pandas as pd
from sklearn.preprocessing import (
    MinMaxScaler,
    RobustScaler,
    StandardScaler,
    LabelBinarizer,
)
from data_prep_tools import *

# 1. Defining variables------------------------------------------------------------------------

# Setting the file path for the predictor variable (Fleshman yeast display expression measure)
predictor_data_path = (
    "antibodyproduction/data/raw_data/fleishman_antibody_expression_data.csv"
)

# Setting the file path for the features (DE-STRESS metrics)
features_data_path = "antibodyproduction/data/raw_data/fleishman_destress_data.csv"

# Setting the file path for the PDB features (DE-STRESS metrics)
pdb_features_data_path = "antibodyproduction/data/raw_data/pdb_destress_data.csv"

# Setting the data output path
processed_data_path = "antibodyproduction/data/processed_data/"

# Setting the data exploration output path
data_exploration_path = "antibodyproduction/analysis/data_exploration/"

# Defining a list of DE-STRESS metrics which are energy field metrics
energy_field_list = [
    "hydrophobic_fitness",
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
    "aggrescan3d_total_value",
    "rosetta_pro_close",
    # "rosetta_omega",
    # "rosetta_total",
    # "rosetta_fa_rep",
    # "evoef2_total",
    # "evoef2_interS_total",
]

# Defining the scaling methods list
scaling_method_list = ["standard", "robust", "minmax"]

# Defining a composition metrics included flag
comp_flag_list = ["comp", "no_comp"]

# Defining a threshold for the spearman correlation coeffient
# in order to remove highly correlated variables
corr_coeff_threshold = 0.7

# Defining a threshold to remove features that have pretty much the same value
constant_features_threshold = 0.8

# Setting variables for train_test_split()
test_size = 0.25
random_state = 42

# 2. Reading in data sets-------------------------------------------------------------------------------

# Reading in Fleishman expression data
raw_predictor_data = pd.read_csv(predictor_data_path)

# Reading in DE-STRESS data
raw_features_data = pd.read_csv(features_data_path)

# Reading in PDB DE-STRESS data
raw_pdb_features_data = pd.read_csv(pdb_features_data_path)

# 3. Processing predictor data (Fleishman expression data)-----------------------------------------------

predictor_data = raw_predictor_data

# Removing blank spaces in the sequences
predictor_data["sequence"] = predictor_data["sequence"].str.replace(" ", "")

# Adding a new field which extracts the design cycle from the design name
# (excluding 4m5.3 as this is the original template antibody)
predictor_data["design_cycle"] = np.where(
    predictor_data["design_name"] == "4m53",
    "",
    predictor_data["design_name"].str.slice(start=0, stop=1),
)

# Adding a new field which extracts the sequence length
predictor_data["seq_len"] = predictor_data["sequence"].str.len()

# Adding a new field which extracts the target antigen
predictor_data["target_antigen"] = np.where(
    predictor_data["design_name"].str.contains("ins"),
    "Insulin",
    np.where(predictor_data["design_name"].str.contains("acp"), "Tuberculosis ACP", ""),
)

# Adding a new field to create an expression bin
predictor_data["expression_bin"] = pd.cut(predictor_data["expression"], 3).astype("str")

# Adding a new field to create an sequence length bin
predictor_data["expression_bin_label"] = np.select(
    [
        predictor_data["expression_bin"].eq("(0.0177, 27.533]"),
        predictor_data["expression_bin"].eq("(27.533, 54.967]"),
        predictor_data["expression_bin"].eq("(54.967, 82.4]"),
    ],
    ["Low", "Medium", "High"],
    default="Unknown",
)

# Adding a new field to create an sequence length bin
predictor_data["length_bin"] = np.select(
    [
        predictor_data["seq_len"].ge(235) & predictor_data["seq_len"].le(240),
        predictor_data["seq_len"].gt(240) & predictor_data["seq_len"].le(245),
        predictor_data["seq_len"].gt(245) & predictor_data["seq_len"].le(250),
        predictor_data["seq_len"].gt(250) & predictor_data["seq_len"].le(255),
        predictor_data["seq_len"].gt(255) & predictor_data["seq_len"].le(260),
        predictor_data["seq_len"].gt(260) & predictor_data["seq_len"].le(265),
    ],
    ["235-240", "240-245", "245-250", "250-255", "255-260", "260-265"],
    default="Unknown",
)

# Removing 4m53 design
predictor_data = predictor_data[predictor_data["design_name"] != "4m53"].reset_index(
    drop=True
)

# Sorting values by design name
processed_predictor_data = predictor_data.sort_values("design_name").reset_index(
    drop=True
)

# Outputting expression data
processed_predictor_data.to_csv(
    processed_data_path + "processed_expression_data.csv", index=False
)


# 4. Processing features (DE-STRESS metrics)----------------------------------------------------------

destress_cols = raw_features_data.columns.to_list()

raw_features_data["top_rank_flag"] = np.where(
    raw_features_data["design_name"].str.contains("rank_001"), 1, 0
)

raw_features_data["relaxed_flag"] = np.where(
    raw_features_data["design_name"].str.contains("unrelaxed"), 0, 1
)

raw_features_data = raw_features_data.loc[
    (
        (raw_features_data["top_rank_flag"] == 1)
        & (raw_features_data["relaxed_flag"] == 0)
    )
].reset_index(drop=True)

raw_features_data.drop(["top_rank_flag", "relaxed_flag"], axis=1, inplace=True)

# Processing Fleishman DE-STRESS data
features_data, design_name_df = process_destress_data(
    raw_features_data,
    energy_field_list=energy_field_list,
    drop_cols=drop_cols,
    data_exploration_path=data_exploration_path,
    mean_std_id="destress_data_num",
)


raw_features_data.to_csv(
    processed_data_path + "fleishman_scfvs_destress_proteome_check.csv", index=False
)

raw_pdb_features_data = raw_pdb_features_data[destress_cols]

raw_pdb_features_data = raw_pdb_features_data.dropna()

# raw_pdb_features_data["top_rank_flag"] = np.where(
#     raw_pdb_features_data["design_name"].str.contains("rank_1"), 1, 0
# )

# raw_pdb_features_data["relaxed_flag"] = np.where(
#     raw_pdb_features_data["design_name"].str.contains("unrelaxed"), 0, 1
# )

# raw_pdb_features_data = raw_pdb_features_data.loc[
#     (
#         (raw_pdb_features_data["top_rank_flag"] == 1)
#         & (raw_pdb_features_data["relaxed_flag"] == 0)
#     )
# ].reset_index(drop=True)

# raw_pdb_features_data.drop(["top_rank_flag", "relaxed_flag"], axis=1, inplace=True)

# Processing PDB DE-STRESS data
pdb_features_data, pdb_design_name_df = process_destress_data(
    raw_pdb_features_data,
    energy_field_list=energy_field_list,
    drop_cols=drop_cols,
    data_exploration_path=data_exploration_path,
    mean_std_id="pdb_destress_data_num",
)

# Removing constant features from Fleishman DE-STRESS features
(
    features_data_constant_removed,
    constant_features,
) = remove_constant_features(
    data=features_data,
    constant_features_threshold=constant_features_threshold,
    output_path=data_exploration_path,
)

# Removing these features from PDB features as well
pdb_features_data_constant_removed = pdb_features_data.drop(
    constant_features, axis=1, inplace=False
)

print("Features dropped because they're constant")
print(constant_features)


# plot_hists_all_columns(
#     data=features_data_constant_removed,
#     column_list=features_data_constant_removed.columns.to_list(),
#     output_path=data_exploration_path,
#     file_name="/pre_scaling_hist_",
# )

# plot_hists_all_columns(
#     data=pdb_features_data_constant_removed,
#     column_list=pdb_features_data_constant_removed.columns.to_list(),
#     output_path=data_exploration_path,
#     file_name="/pdb_pre_scaling_hist_",
# )


# Adding design name back on for the train and test split
features_data_constant_removed = pd.concat(
    [features_data_constant_removed, design_name_df], axis=1
)
pdb_features_data_constant_removed = pd.concat(
    [pdb_features_data_constant_removed, pdb_design_name_df], axis=1
)

# 5. Splitting data up into train and test and scaling----------------------------------------------------------

# Scaling features
for scaling_method in scaling_method_list:
    for comp_flag in comp_flag_list:
        scaled_processed_data_path = (
            processed_data_path + scaling_method + "/" + comp_flag + "/"
        )
        scaled_data_exporation_path = (
            data_exploration_path + scaling_method + "/" + comp_flag + "/"
        )

        # Train and test splits
        X_train, X_test, y_train, y_test = train_test_split(
            features_data_constant_removed,
            processed_predictor_data,
            test_size=test_size,
            random_state=random_state,
            stratify=processed_predictor_data["expression_bin"],
        )

        # Resetting index for PDB data
        pdb_features_data_processed = pdb_features_data_constant_removed.reset_index(
            drop=True
        )

        # Removing composition metrics so we have training sets without them
        if comp_flag == "no_comp":
            X_train = X_train[
                X_train.columns.drop(list(X_train.filter(regex="composition")))
            ]
            X_test = X_test[
                X_test.columns.drop(list(X_test.filter(regex="composition")))
            ]

            pdb_features_data_processed = pdb_features_data_processed[
                pdb_features_data_processed.columns.drop(
                    list(pdb_features_data_processed.filter(regex="composition"))
                )
            ]

        # Sorting all the data sets by design name
        X_train = X_train.sort_values("design_name")
        X_test = X_test.sort_values("design_name")
        y_train = y_train.sort_values("design_name")
        y_test = y_test.sort_values("design_name")
        pdb_features_data_processed.sort_values("design_name")

        # Binarizing the y labels
        label_binarizer = LabelBinarizer().fit(
            y_train["expression_bin_label"].to_numpy()
        )

        y_train_binary = label_binarizer.transform(
            y_train["expression_bin_label"].to_numpy()
        )
        with open(scaled_processed_data_path + "y_train_binary.npy", "wb") as f:
            np.save(f, y_train_binary)

        y_test_binary = label_binarizer.transform(
            y_test["expression_bin_label"].to_numpy()
        )

        with open(scaled_processed_data_path + "y_test_binary.npy", "wb") as f:
            np.save(f, y_test_binary)

        # Exporting data sets as csv files
        X_train.to_csv(scaled_processed_data_path + "X_train.csv", index=False)
        X_test.to_csv(scaled_processed_data_path + "X_test.csv", index=False)
        y_train.to_csv(scaled_processed_data_path + "y_train.csv", index=False)
        y_test.to_csv(scaled_processed_data_path + "y_test.csv", index=False)

        # Saving design name order for training data
        X_train["design_name"].to_csv(
            scaled_processed_data_path + "train_design_name_order.csv", index=False
        )

        # Saving design name order for training data
        X_test["design_name"].to_csv(
            scaled_processed_data_path + "test_design_name_order.csv", index=False
        )

        # Saving design name order for PDB data
        pdb_features_data_processed["design_name"].to_csv(
            scaled_processed_data_path + "pdb_design_name_order.csv", index=False
        )

        # Now removing design name from X_train and X_test
        X_train.drop(
            [
                "design_name",
            ],
            axis=1,
            inplace=True,
        )

        # Dropping design name
        X_test.drop(
            [
                "design_name",
            ],
            axis=1,
            inplace=True,
        )

        # Dropping design name
        pdb_features_data_processed.drop(
            [
                "design_name",
            ],
            axis=1,
            inplace=True,
        )

        # Scaling data
        if scaling_method == "minmax":
            # Scaling the data with min max scaler
            scaler = MinMaxScaler().fit(X_train)

        elif scaling_method == "standard":
            # Scaling the data with standard scaler scaler
            scaler = StandardScaler().fit(X_train)

        elif scaling_method == "robust":
            # Scaling the data with robust scaler
            scaler = RobustScaler().fit(X_train)

        # Transforming the data
        X_train_scaled = pd.DataFrame(
            scaler.transform(X_train),
            columns=X_train.columns,
        )
        # Transforming test data with same scaler
        X_test_scaled = pd.DataFrame(
            scaler.transform(X_test),
            columns=X_train.columns,
        )

        # Transforming PDB data with same scaler
        pdb_scaled = pd.DataFrame(
            scaler.transform(pdb_features_data_processed),
            columns=X_train.columns,
        )

        # Calculating mean and std of features
        features_mean_std(
            data=X_train_scaled,
            output_path=scaled_data_exporation_path,
            id="destress_data_scaled",
        )

        (
            X_train_drop_high_corr,
            drop_cols_high_corr,
        ) = remove_highest_correlators(
            data=X_train_scaled,
            corr_coeff_threshold=corr_coeff_threshold,
            output_path=scaled_data_exporation_path,
        )

        print("Features dropped because of high correlation")
        print(drop_cols_high_corr)

        # Remove these features from test as well
        X_test_drop_high_corr = X_test_scaled.drop(
            drop_cols_high_corr,
            axis=1,
            inplace=False,
        )
        pdb_drop_high_corr = pdb_scaled.drop(
            drop_cols_high_corr,
            axis=1,
            inplace=False,
        )

        # plot_hists_all_columns(
        #     data=X_train_drop_high_corr,
        #     column_list=X_train_drop_high_corr.columns.to_list(),
        #     output_path=scaled_data_exporation_path,
        #     file_name="/post_scaling_hist_",
        # )
        # plot_hists_all_columns(
        #     data=pdb_drop_high_corr,
        #     column_list=pdb_drop_high_corr.columns.to_list(),
        #     output_path=scaled_data_exporation_path,
        #     file_name="/pdb_post_scaling_hist_",
        # )

        X_train_drop_high_corr.to_csv(
            scaled_processed_data_path + "X_train_scaled.csv",
            index=False,
        )

        X_test_drop_high_corr.to_csv(
            scaled_processed_data_path + "X_test_scaled.csv",
            index=False,
        )

        pdb_drop_high_corr.to_csv(
            scaled_processed_data_path + "pdb_scaled.csv",
            index=False,
        )
