# This script provides helper functions which are used in the data prep scripts

# 0. Loading in packages and defining custom functions--------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.preprocessing import MinMaxScaler, StandardScaler, RobustScaler
import pickle

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


# Defining a fucntion add destress summary columns
def adding_destress_summary_cols(destress_data):

    # Adding a new field to create a dssp bin
    destress_data["dssp_bin"] = np.select(
        [
            destress_data["ss_prop_alpha_helix"].gt(0.5),
            destress_data["ss_prop_beta_bridge"].gt(0.5),
            destress_data["ss_prop_beta_strand"].gt(0.5),
            destress_data["ss_prop_3_10_helix"].gt(0.5),
            destress_data["ss_prop_pi_helix"].gt(0.5),
            destress_data["ss_prop_hbonded_turn"].gt(0.5),
            destress_data["ss_prop_bend"].gt(0.5),
            destress_data["ss_prop_loop"].gt(0.5),
        ],
        [
            "Alpha Helix",
            "Beta Bridge",
            "Beta Strand",
            "3 10 Helix",
            "Pi Helix",
            "Hbond Turn",
            "Bend",
            "Loop",
        ],
        default="Mixed",
    )

    # Adding a new field to create a isoelectric point bin
    destress_data["isoelectric_point_bin"] = np.select(
        [
            destress_data["isoelectric_point"].lt(6),
            destress_data["isoelectric_point"].ge(6)
            & destress_data["isoelectric_point"].le(8),
            destress_data["isoelectric_point"].gt(8),
        ],
        [
            "Less than 6",
            "Between 6 and 8",
            "Greater than 8",
        ],
        default="Unknown",
    )

    # Adding a new field to create a packing density bin
    destress_data["packing_density_bin"] = np.select(
        [
            destress_data["packing_density"].lt(40),
            destress_data["packing_density"].ge(40)
            & destress_data["packing_density"].lt(60),
            destress_data["packing_density"].ge(60)
            & destress_data["packing_density"].lt(80),
            destress_data["packing_density"].ge(80),
        ],
        [
            "Less than 40",
            "Between 40 and 60",
            "Between 60 and 80",
            "Greater than 80",
        ],
        default="Unknown",
    )

    # Adding a new field to create a packing density bin
    destress_data["aggrescan3d_avg_bin"] = np.select(
        [
            destress_data["aggrescan3d_avg_value"].lt(-2),
            destress_data["aggrescan3d_avg_value"].ge(-2)
            & destress_data["aggrescan3d_avg_value"].lt(0),
            destress_data["aggrescan3d_avg_value"].ge(0)
            & destress_data["aggrescan3d_avg_value"].lt(2),
            destress_data["aggrescan3d_avg_value"].ge(2),
        ],
        [
            "Less than -2",
            "Between -2 and 0",
            "Between 0 and 2",
            "Greater than 2",
        ],
        default="Unknown",
    )

    return destress_data


# Adding a function to join on the af2db uniprot data and make some organism group columns
def adding_af2db_uniprot_columns(
    destress_data,
    af2db_uniprot_data,
    af2db_cluster_data,
    organism_group_dict,
):

    # Extracting the uniprot id from the design name in the DE-STRESS data so that we can join on the AF2DB and uniprot data
    destress_data["uniprot_id"] = destress_data["design_name"].str.split("-").str[1]

    # Joining on af2db and uniprot data
    destress_af2db_uniprot_data = pd.merge(
        destress_data,
        af2db_uniprot_data[
            [
                "uniprot_id",
                "organism_scientific_name",
                "subcellular_location",
                "uniprot_description",
            ]
        ],
        on="uniprot_id",
        how="left",
    )

    print(destress_af2db_uniprot_data)

    # Joining on af2db cluster data
    destress_af2db_uniprot_data = pd.merge(
        destress_af2db_uniprot_data,
        af2db_cluster_data[
            [
                "uniprot_id",
                "cluster_representative",
            ]
        ],
        on="uniprot_id",
        how="left",
    )

    print(destress_af2db_uniprot_data)

    # Adding a new field to create an organism group
    destress_af2db_uniprot_data["organism_group"] = np.select(
        [
            destress_af2db_uniprot_data["organism_scientific_name"].isin(
                organism_group_dict["animal"]
            ),
            destress_af2db_uniprot_data["organism_scientific_name"].isin(
                organism_group_dict["fungi"]
            ),
            destress_af2db_uniprot_data["organism_scientific_name"].isin(
                organism_group_dict["bacteria"]
            ),
            destress_af2db_uniprot_data["organism_scientific_name"].isin(
                organism_group_dict["plant"]
            ),
            destress_af2db_uniprot_data["organism_scientific_name"].isin(
                organism_group_dict["protozoan"]
            ),
            destress_af2db_uniprot_data["organism_scientific_name"].isin(
                organism_group_dict["archaea"]
            ),
            destress_af2db_uniprot_data["organism_scientific_name"].isin(
                organism_group_dict["other"]
            ),
        ],
        [
            "Animal",
            "Fungi",
            "Bacteria",
            "Plant",
            "Protozoan",
            "Archaea",
            "Other",
        ],
        default="Unknown",
    )

    print(destress_af2db_uniprot_data.value_counts("organism_group"))

    # Adding a new field to create an organism group
    destress_af2db_uniprot_data["organism_group2"] = np.select(
        [
            destress_af2db_uniprot_data["organism_scientific_name"].isin(
                organism_group_dict["animal"]
                + organism_group_dict["protozoan"]
                + organism_group_dict["fungi"]
                + organism_group_dict["plant"]
            ),
            destress_af2db_uniprot_data["organism_scientific_name"].isin(
                organism_group_dict["bacteria"] + organism_group_dict["archaea"]
            ),
            destress_af2db_uniprot_data["organism_scientific_name"].isin(
                organism_group_dict["other"]
            ),
        ],
        [
            "Eukaryotes",
            "Prokaryotes",
            "Other",
        ],
        default="Unknown",
    )

    print(destress_af2db_uniprot_data.value_counts("organism_group2"))

    return destress_af2db_uniprot_data


# Defining a function to remove destress features that have low variance, are non numeric or defined in drop_columns_list
def filtering_destress_metrics(
    destress_data, drop_cols_list, data_exploration_path, constant_features_threshold
):

    destress_columns_full = destress_data.columns.to_list()

    # Dropping columns that have been defined manually
    destress_data = destress_data.drop(drop_cols_list, axis=1)

    # Dropping columns that are not numeric
    destress_data_num = destress_data.select_dtypes([np.number]).reset_index(drop=True)

    data_corr = stats.spearmanr(destress_data_num)
    data_corr_df = pd.DataFrame(
        data_corr[0],
        columns=destress_data_num.columns.to_list(),
        index=destress_data_num.columns.to_list(),
    )

    data_corr_df.to_csv(data_exploration_path + "corr.csv")

    # Dropping composition metrics
    destress_data_num = destress_data_num[
        destress_data_num.columns.drop(
            list(destress_data_num.filter(regex="composition"))
        )
    ]

    # Printing columns that are dropped because they are not numeric
    destress_columns_num = destress_data_num.columns.to_list()
    dropped_cols_non_num = set(destress_columns_full) - set(destress_columns_num)

    print("Features dropped because they're not numeric or in drop_cols_list")
    print(dropped_cols_non_num)

    # Calculating mean and std of features
    features_mean_std(
        data=destress_data_num,
        output_path=data_exploration_path,
        id="destress_data_num",
    )

    (
        destress_data_filtered,
        constant_features,
    ) = remove_constant_features(
        data=destress_data_num,
        constant_features_threshold=constant_features_threshold,
        output_path=data_exploration_path,
    )

    print("Features dropped because they're constant")
    print(constant_features)

    return destress_data_filtered


# Defining a function to scale the data with different scaling methods and remove the highest correlated features
def scale_destress_data_remove_high_corr(
    destress_data,
    scaling_method,
    data_exploration_path,
    data_output_path,
    corr_coeff_threshold,
):

    data_exploration_scaled_output_path = data_exploration_path + scaling_method + "/"
    data_scaled_output_path = data_output_path + scaling_method + "/"

    if scaling_method == "minmax":
        # Scaling the data with min max scaler
        scaler = MinMaxScaler().fit(destress_data)
        # Making a dump of the fitted scaler
        with open(data_scaled_output_path + "minmax_scaler.pkl", "wb") as file:
            pickle.dump(scaler, file)

    elif scaling_method == "standard":
        # Scaling the data with standard scaler scaler
        scaler = StandardScaler().fit(destress_data)
        # Making a dump of the fitted scaler
        with open(data_scaled_output_path + "standard_scaler.pkl", "wb") as file:
            pickle.dump(scaler, file)

    elif scaling_method == "robust":
        # Scaling the data with robust scaler
        scaler = RobustScaler().fit(destress_data)
        # Making a dump of the fitted scaler
        with open(data_scaled_output_path + "robust_scaler.pkl", "wb") as file:
            pickle.dump(scaler, file)

    # Transforming the data
    destress_data_scaled = pd.DataFrame(
        scaler.transform(destress_data),
        columns=destress_data.columns,
    )

    print("Scaled data column order")
    print(destress_data_scaled.columns.to_list())

    # Calculating mean and std of features
    features_mean_std(
        data=destress_data_scaled,
        output_path=data_exploration_scaled_output_path,
        id="destress_data_scaled",
    )

    (
        destress_data_remove_high_corr,
        drop_cols_high_corr,
    ) = remove_highest_correlators(
        data=destress_data_scaled,
        corr_coeff_threshold=corr_coeff_threshold,
        output_path=data_exploration_scaled_output_path,
    )

    print("Features dropped because of high correlation")
    print(drop_cols_high_corr)

    # plot_hists_all_columns(
    #     data=destress_data_remove_high_corr,
    #     column_list=destress_data_remove_high_corr.columns.to_list(),
    #     output_path=data_exploration_scaled_output_path,
    #     file_name="/post_scaling_hist_",
    # )

    destress_data_remove_high_corr.to_csv(
        data_scaled_output_path + "processed_destress_data_scaled.csv",
        index=False,
    )

    return destress_data_remove_high_corr


# Defining a function to process the af2 data
def process_af2_data(
    raw_destress_data,
    af2db_uniprot_data,
    af2db_cluster_data,
    data_exploration_path,
    data_output_path,
    missing_val_threshold,
    organism_group_dict,
    energy_field_list,
    labels,
    drop_cols_list,
    constant_features_threshold,
    scaling_method_list,
    corr_coeff_threshold,
):

    # Removing features that have missing value prop greater than threshold
    destress_data, dropped_cols_miss_vals = remove_missing_val_features(
        data=raw_destress_data,
        output_path=data_exploration_path,
        threshold=missing_val_threshold,
    )

    print("Columns dropped because of missing values")
    print(dropped_cols_miss_vals)

    #  Calculating total number of structures that DE-STRESS ran for
    num_structures = destress_data.shape[0]

    # Now removing any rows that have missing values
    destress_data = destress_data.dropna(axis=0).reset_index(drop=True)

    # Calculating number of structures in the data set after removing missing values
    num_structures_missing_removed = destress_data.shape[0]

    # Calculating how many structures are left after removing those with missing values for the DE-STRESS metrics.
    print(
        "DE-STRESS ran for "
        + str(num_structures)
        + " AF2 structures in total and after removing missing values there are "
        + str(num_structures_missing_removed)
        + " structures remaining in the data set. This means "
        + str(100 * (round((num_structures_missing_removed / num_structures), 4)))
        + "% of the protein structures are covered in this data set."
    )

    # Adding DE-STRESS summary columns
    destress_data = adding_destress_summary_cols(destress_data=destress_data)

    # Joining on af2b and uniprot data and making organism columns
    destress_af2db_uniprot_data = adding_af2db_uniprot_columns(
        destress_data=destress_data,
        af2db_uniprot_data=af2db_uniprot_data,
        af2db_cluster_data=af2db_cluster_data,
        organism_group_dict=organism_group_dict,
    )

    # Normalising energy field values by the number of residues
    destress_af2db_uniprot_data.loc[
        :,
        energy_field_list,
    ] = destress_af2db_uniprot_data.loc[
        :,
        energy_field_list,
    ].div(
        destress_af2db_uniprot_data["num_residues"], axis=0
    )

    # Saving labels
    labels_df = save_destress_labels(
        data=destress_af2db_uniprot_data,
        labels=labels,
        output_path=data_output_path,
        file_path="labels",
    )
    print(labels_df)

    # Filtering destress columns from drop cols list, low variance/constant and non numeric
    destress_data_filtered = filtering_destress_metrics(
        destress_data=destress_af2db_uniprot_data,
        drop_cols_list=drop_cols_list,
        data_exploration_path=data_exploration_path,
        constant_features_threshold=constant_features_threshold,
    )

    # Scaling the data sets
    for scaling_method in scaling_method_list:
        print(scaling_method)
        scaled_destress_data = scale_destress_data_remove_high_corr(
            destress_data=destress_data_filtered,
            scaling_method=scaling_method,
            data_exploration_path=data_exploration_path,
            data_output_path=data_output_path,
            corr_coeff_threshold=corr_coeff_threshold,
        )

        print(scaled_destress_data)
