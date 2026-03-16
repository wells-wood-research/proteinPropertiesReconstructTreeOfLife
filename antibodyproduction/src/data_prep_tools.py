# This script defines helper functions which are used throughout the
# data prep script.

# 0. Importing packages----------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
import pickle


# 1. Defining helper functions---------------------------------------------------


# Defining a function to process the DE-STRESS data
def process_destress_data(
    data, energy_field_list, drop_cols, data_exploration_path, mean_std_id
):
    features_data = data

    # Extracting the design name
    features_data["design_name"] = features_data["design_name"].str.split("_").str[0]

    # Removing 4m53 design
    features_data = features_data[features_data["design_name"] != "4m53"]

    # Sorting the values by design name
    features_data = features_data.sort_values("design_name").reset_index(drop=True)

    # Extracting design name so we can join it back on later
    design_name_df = features_data["design_name"]

    # Normalising energy field values by the number of residues
    features_data.loc[
        :,
        energy_field_list,
    ] = features_data.loc[
        :,
        energy_field_list,
    ].div(features_data["num_residues"], axis=0)

    features_data_columns = features_data.columns.to_list()

    # Dropping columns that have been defined manually
    features_data_filt = features_data.drop(drop_cols, axis=1)

    # Dropping columns that are not numeric
    features_data_num = features_data_filt.select_dtypes([np.number]).reset_index(
        drop=True
    )

    # Printing columns that are dropped because they are not numeric
    destress_columns_num = features_data_num.columns.to_list()
    dropped_cols_non_num = set(features_data_columns) - set(destress_columns_num)

    print("Dropping non numeric columns and columns specified in drop cols")
    print(dropped_cols_non_num)

    # Calculating mean and std of features
    features_mean_std(
        data=features_data_num,
        output_path=data_exploration_path,
        id=mean_std_id,
    )

    return features_data_num, design_name_df


# Defining a function to compute mean and std of features
def features_mean_std(data, output_path, id):
    data_std = data.std().sort_values(ascending=False)
    data_mean = data.mean().sort_values(ascending=False)

    data_std.to_csv(output_path + "data_std_" + id + ".csv")
    data_mean.to_csv(output_path + "data_mean_" + id + ".csv")


# Defining a function to plot histograms for all columns in a data set
def plot_hists_all_columns(data, column_list, output_path, file_name):
    for col in column_list:
        plt.hist(data=data, x=col, bins=50)
        # plt.hist(data=data, x=col, histtype="step")
        plt.savefig(
            output_path + file_name + col + ".png",
            bbox_inches="tight",
            dpi=300,
        )
        plt.close()


def remove_highest_correlators(data, corr_coeff_threshold, output_path):
    drop_cols_high_corr = []

    data_corr = stats.spearmanr(data)
    data_corr_df = pd.DataFrame(
        data_corr[0], columns=data.columns.to_list(), index=data.columns.to_list()
    )
    data_corr_df.to_csv(output_path + "corr_matrix_before.csv", index=False)
    data_corr_df_abs = data_corr_df.abs()
    data_cutoff_count = (
        data_corr_df_abs[data_corr_df_abs > corr_coeff_threshold]
        .count()
        .sort_values(ascending=False)
    )

    while data_cutoff_count.max() > 1:
        data_corr = stats.spearmanr(data)
        data_corr_df = pd.DataFrame(
            data_corr[0], columns=data.columns.to_list(), index=data.columns.to_list()
        )
        data_corr_df_abs = data_corr_df.abs()
        data_cutoff_count = (
            data_corr_df_abs[data_corr_df_abs > corr_coeff_threshold]
            .count()
            .sort_values(ascending=False)
        )
        data_cutoff_count_max = data_cutoff_count.index.values[0]
        drop_cols_high_corr.append(data_cutoff_count_max)

        data = data.drop(data_cutoff_count_max, axis=1)

    data_corr_df.to_csv(output_path + "corr_matrix_after.csv", index=False)

    data_new = data

    return data_new, drop_cols_high_corr


# Defining a function to filter out further features
def remove_constant_features(data, constant_features_threshold, output_path):
    prop_max_value_count_list = []

    for col in data.columns.to_list():
        prop_max_value_count = np.max(round(data[col], 2).value_counts(normalize=True))

        prop_max_value_count_list.append(prop_max_value_count)

    prop_max_value_count_df = pd.DataFrame(
        dict(
            zip(
                ["features", "prop_max_value_count"],
                [data.columns.to_list(), prop_max_value_count_list],
            )
        )
    )

    prop_max_value_count_df.sort_values(
        by="prop_max_value_count", inplace=True, ascending=False
    )

    prop_max_value_count_df.to_csv(
        output_path + "prop_max_value_count_df.csv", index=False
    )

    constant_features = prop_max_value_count_df["features"][
        prop_max_value_count_df["prop_max_value_count"] > constant_features_threshold
    ].to_list()

    data_new = data.drop(constant_features, axis=1).reset_index(drop=True)

    return data_new, constant_features


# Function to create train and test split
def create_train_test_data(
    predictor_data: pd.DataFrame,
    features_data: pd.DataFrame,
    test_size: float,
    random_state: int,
    train_output_path: str,
    test_output_path: str,
) -> None:
    # Getting rid of 4m5.3 for now
    features_data = features_data[features_data["design_name"] != "4m53"]
    predictor_data = predictor_data[predictor_data["design_name"] != "4m53"]

    # Train and test splits--------------------------------------
    X_train, X_test, y_train, y_test = train_test_split(
        features_data,
        predictor_data,
        test_size=test_size,
        random_state=random_state,
        stratify=predictor_data["expression_bin"],
    )

    # Sorting all the data sets by design name
    X_train = X_train.sort_values("design_name")
    X_test = X_test.sort_values("design_name")
    y_train = y_train.sort_values("design_name")
    y_test = y_test.sort_values("design_name")

    # Exporting data sets as csv files
    X_train.to_csv(train_output_path + "X_train.csv", index=False)
    X_test.to_csv(test_output_path + "X_test.csv", index=False)
    y_train.to_csv(train_output_path + "y_train.csv", index=False)
    y_test.to_csv(test_output_path + "y_test.csv", index=False)

    # Saving design name order for training data
    X_train["design_name"].to_csv(
        train_output_path + "train_data_design_name_order.csv", index=False
    )

    # Saving design name order for training data
    X_test["design_name"].to_csv(
        test_output_path + "test_data_design_name_order.csv", index=False
    )

    # Now removing design name from X_train and X_test
    X_train.drop(
        [
            "design_name",
        ],
        axis=1,
        inplace=True,
    )

    X_test.drop(
        [
            "design_name",
        ],
        axis=1,
        inplace=True,
    )

    # Scaling X_train and X_test in exactly the same way

    # Fitting the scalers
    X_scaler = MinMaxScaler().fit(X_train)

    pickle.dump(
        X_scaler,
        open(
            "models/scalers/trained_min_max_scaler.p",
            "wb",
        ),
    )

    # Transforming the data
    X_train_scaled = pd.DataFrame(X_scaler.transform(X_train), columns=X_train.columns)
    X_test_scaled = pd.DataFrame(X_scaler.transform(X_test), columns=X_test.columns)

    # Removing composition metrics so we have training sets without them
    X_train_scaled_no_comp = X_train_scaled[
        X_train_scaled.columns.drop(list(X_train_scaled.filter(regex="composition")))
    ]
    X_test_scaled_no_comp = X_test_scaled[
        X_test_scaled.columns.drop(list(X_test_scaled.filter(regex="composition")))
    ]

    # Exporting data sets as csv files
    X_train_scaled.to_csv(train_output_path + "X_train_scaled.csv", index=False)
    X_train_scaled_no_comp.to_csv(
        train_output_path + "X_train_scaled_no_comp.csv", index=False
    )
    X_test_scaled.to_csv(test_output_path + "X_test_scaled.csv", index=False)
    X_test_scaled_no_comp.to_csv(
        test_output_path + "X_test_scaled_no_comp.csv", index=False
    )
