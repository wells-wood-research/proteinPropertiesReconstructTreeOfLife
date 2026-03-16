# This script provides helper functions for feature selection.

# 0. Importing packages-------------------------------------------
import numpy as np
import pandas as pd
from sklearn.feature_selection import mutual_info_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel


# 1. Defining helper functions------------------------------------
# Defining a function to select features by using mutual infiormation
# of each feature against the expression output variable
def feature_select_mi(X_train, y_train_pred_var):
    # Using mutual information for feature selection
    feature_mutual_info = mutual_info_classif(
        X_train, y_train_pred_var, random_state=42
    )

    mutual_info_dict = {
        "feature_names": X_train.columns,
        "mutual_info": feature_mutual_info,
    }

    # Mutual Info data frame
    feature_mutual_info_df = pd.DataFrame(mutual_info_dict)

    # Filtering for features with mutual information greater than or equt to average
    selected_features = feature_mutual_info_df[
        feature_mutual_info_df["mutual_info"].ge(
            feature_mutual_info_df["mutual_info"].mean()
        )
    ].reset_index(drop=True)

    selected_features.sort_values(by="mutual_info", inplace=True, ascending=False)

    return selected_features


# Defining a function to feature select using a random forrest
def feature_select_rf(X_train, y_train_pred_var):
    sel = SelectFromModel(
        RandomForestClassifier(
            n_estimators=1000, class_weight="balanced", random_state=42
        )
    )
    sel.fit(X_train, y_train_pred_var)
    selected_features = X_train.columns[(sel.get_support())]
    feature_importance = sel.estimator_.feature_importances_[(sel.get_support())]
    selected_features_dict = {
        "feature_names": selected_features,
        "feature_importance": feature_importance,
    }

    selected_features = pd.DataFrame(selected_features_dict)

    selected_features.sort_values(
        by="feature_importance", inplace=True, ascending=False
    )

    return selected_features
