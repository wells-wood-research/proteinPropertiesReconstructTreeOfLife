# This script performs dimensionality reduction analysis on the
# Fleishman antibody expression data.

# 0. Importing packages-----------------------------------------------------------
import numpy as np
import pandas as pd
from dim_red_tools import *

# 1. Defining variables------------------------------------------------------------

# Defining the scaling methods list
scaling_method_list = ["standard", "robust", "minmax"]
# scaling_method_list = ["robust"]

# Defining a composition metrics included flag
comp_flag_list = ["comp", "no_comp"]
# comp_flag_list = ["comp"]

# Defining a list of feature selection methods
feature_selection_list = ["mi", "rf"]
# feature_selection_list = ["mi"]

# Defining feature selection path
feature_selection_path = "antibodyproduction/feature_selection/"

# Defining pca output path
pca_output_path = "antibodyproduction/analysis/dim_red/pca/"

# Defining the number of components for PCA
pca_num_components = 10

# Defining list of dim ids
dim_ids_list = []
for i in range(0, pca_num_components):
    dim_ids_list.append("dim" + str(i))

# Defining hover data for plotly
hover_data = ["design_name", "dim0", "dim1"]

# Creating a color palette
# palette = sns.color_palette(["#0173b2", "#d55e00", "#029e73", "#cc78bc"], 4)
palette = sns.color_palette(["#17becf", "#e6a514", "#8c564b", "#3f2b96"], 4)

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
        pca_scaler_output_path = (
            pca_output_path + scaling_method + "/" + comp_flag + "/"
        )
        pdb_design_name_ids_path = train_data_path + "pdb_design_name_order.csv"

        # 2. Reading in data------------------------------------------------------------

        y_train = pd.read_csv(train_data_path + "y_train.csv")
        y_train_pred_var = y_train["expression_bin_label"]
        pdb_design_name_ids = pd.read_csv(pdb_design_name_ids_path)
        pdb_expression_label_df = pdb_design_name_ids
        pdb_expression_label_df["expression_bin_label"] = "SAbDab"

        print(pdb_expression_label_df)

        for feature_selection in feature_selection_list:
            X_train_scaled = pd.read_csv(train_data_path + "X_train_scaled.csv")
            pdb_scaled = pd.read_csv(train_data_path + "pdb_scaled.csv")

            selected_features = pd.read_csv(
                feature_selection_scaler_path
                + "selected_features_"
                + feature_selection
                + ".csv"
            )
            X_train_scaled_filt = X_train_scaled[selected_features["feature_names"]]
            pdb_scaled_filt = pdb_scaled[selected_features["feature_names"]]

            # 3. Performing PCA--------------------------------------------------------------------------

            # pca_var_explained(
            #     data=X_train_scaled_filt,
            #     n_components=pca_num_components,
            #     file_name="pca_var_explained_" + feature_selection,
            #     output_path=pca_scaler_output_path,
            # )

            # pca_transformed_data, pca_model = perform_pca(
            #     data=X_train_scaled_filt,
            #     labels_df=y_train[["design_name", "expression_bin_label"]],
            #     n_components=pca_num_components,
            #     output_path=pca_scaler_output_path,
            #     file_path="pca_transformed_data_" + feature_selection,
            #     components_file_path="comp_contrib_" + feature_selection + "_",
            # )

            # pdb_pca_transformed_data = pca_model.transform(pdb_scaled_filt)

            # # Converting to data frame and renaming columns
            # pdb_pca_transformed_data = pd.DataFrame(pdb_pca_transformed_data)

            # pdb_pca_transformed_data.columns = dim_ids_list

            # # Adding the labels back
            # pdb_pca_transformed_data = pd.concat(
            #     [pdb_expression_label_df, pdb_pca_transformed_data],
            #     axis=1,
            # )

            # pca_transformed_data_joined = pd.concat(
            #     [pca_transformed_data, pdb_pca_transformed_data], axis=0
            # ).reset_index(drop=True)

            # print(pca_transformed_data_joined)

            joined_X_data = pd.concat(
                [X_train_scaled_filt, pdb_scaled_filt], axis=0
            ).reset_index(drop=True)

            print(joined_X_data)

            joined_y_data = pd.concat(
                [
                    y_train[["design_name", "expression_bin_label"]],
                    pdb_expression_label_df,
                ],
                axis=0,
            ).reset_index(drop=True)

            print(joined_y_data)

            var_explained_df = pca_var_explained(
                data=joined_X_data,
                n_components=pca_num_components,
                file_name="pca_var_explained_" + feature_selection,
                output_path=pca_scaler_output_path,
            )

            pca_transformed_data, pca_model = perform_pca(
                data=joined_X_data,
                labels_df=joined_y_data,
                n_components=pca_num_components,
                output_path=pca_scaler_output_path,
                file_path="pca_transformed_data_" + feature_selection,
                components_file_path="comp_contrib_" + feature_selection + "_",
            )

            # 5. Plotting 2d spaces---------------------------------------------------------------------

            # Setting theme for plots
            sns.set_style("whitegrid")
            # cmap = sns.color_palette("colorblind", as_cmap=True)

            plot_pca_plotly(
                pca_data=pca_transformed_data,
                x="dim0",
                y="dim1",
                z="dim2",
                color="expression_bin_label",
                hover_data=hover_data,
                legend_title="Expression",
                opacity=0.7,
                size=10,
                output_path=pca_scaler_output_path,
                file_name="pca_embedding_expression_" + feature_selection + ".html",
            )

            plot_latent_space_2d(
                data=pca_transformed_data,
                var_explained_data=var_explained_df,
                x="dim0",
                y="dim1",
                axes_prefix="PC",
                legend_title="Expression Level",
                hue="expression_bin_label",
                hue_order=["Low", "Medium", "High", "SAbDab"],
                # style=var,
                alpha=0.7,
                s=150,
                palette=palette,
                output_path=pca_scaler_output_path,
                file_name="pca_embedding_expression_" + feature_selection,
            )

            spectral_plot(
                pca_data=pca_transformed_data,
                group_var="expression_bin_label",
                hue_order=["Low", "Medium", "High", "SAbDab"],
                value_var_list=dim_ids_list,
                title="",
                legend_title="",
                output_path=pca_scaler_output_path,
                file_name="spectral_plot_" + feature_selection,
            )
