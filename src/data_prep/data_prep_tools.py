# This script provides helper functions which are used in the data prep scripts

# 0. Loading in packages and defining custom functions--------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# 1. Removing features-----------------------------------------------------------------------------


# Defining a function which removes features with missing values
# if the proportion is greater than a threshold
def remove_missing_val_features(data, output_path, threshold):
    data_missing_count = data.isnull().sum()
    data_missing_count_df = pd.DataFrame(data_missing_count, columns=["num_missing"])
    data_missing_count_df["prop_missing"] = (
        data_missing_count_df["num_missing"] / data.shape[0]
    )
    data_missing_count_df.to_csv(output_path + "destress_data_missing_count_df.csv")

    features_to_remove = data_missing_count_df[
        data_missing_count_df["prop_missing"] > threshold
    ].index.values

    new_data = data.drop(features_to_remove, axis=1)

    return new_data, features_to_remove


# Defining a function which saves csv files of the labels
def save_destress_labels(data, labels, output_path, file_path):
    data_filt = data[labels]
    data_filt.to_csv(output_path + file_path + ".csv", index=False)

    return data_filt


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
