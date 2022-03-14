#!/usr/bin/env python3
# Author: Norbert Moldován

import os
import sys
import argparse
import pandas as pd

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFECV
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.model_selection import GridSearchCV as GCV
from sklearn.model_selection import LeaveOneOut as LOO

from sklearn.neighbors import KNeighborsRegressor as KNR
from sklearn.linear_model import LogisticRegression as LoR
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier as MLP

from sklearn.metrics import precision_score, recall_score, accuracy_score
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.metrics import roc_auc_score, plot_roc_curve, roc_curve
from sklearn.metrics import make_scorer

from skbio.stats.composition import clr
import numpy as np
import matplotlib
#matplotlib.use('Agg') # Matplotlib can't use interactive backend on HPC.
import matplotlib.pyplot as plt
import seaborn as sns

pd.set_option("display.max_colwidth", None)  # Print full column with for pd.Df.
# Parsing the arguments + some nice help text.We
def argumentparsing():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in",
                        "--input",
                        dest="InputFile",
                        type=str,
                        required=False,
                        help="The path to the imput data. Use this, if "
                        "you want to split the data for validation."
                        "Can't be used together with training/predict inputs!")
    parser.add_argument("-tr",
                        "--training",
                        dest="TrainingFile",
                        type=str,
                        required=False,
                        help="The path to the training data. "
                        "Can't be used together with input!")
    parser.add_argument("-pr",
                        "--predict",
                        dest="PredictFile",
                        type=str,
                        required=False,
                        help="The path to the prediction data. "
                        "Can't be used together with input!")
    parser.add_argument("-o",
                        "--output",
                        dest="OutputFolder",
                        type=str,
                        required=False,
                        help="Output path and prefix.")
    parser.add_argument("-m",
                        "--models",
                        dest="Models",
                        type=str,
                        required=True,
                        help="Models selected for prediction.")
    parser.add_argument("-y",
                        "--target",
                        dest="TargetColumn",
                        type=str,
                        required=True,
                        help="Which column to use as a target? If none is set "
                        "or the set value is not in the files 'Groups' or "
                        "'WhichGroup' will be used.")
    parser.add_argument("-clrt",
                        "--clrt",
                        dest="useClrt",
                        action="store_true",
                        required=False,
                        help="Use the central log transform to Transform "
                        "the dataset.")
    parser.add_argument("-ign",
                        "--ignore_col",
                        dest="ignore_col",
                        type=str,
                        required=False,
                        help="Ignore column.")
    parser.add_argument("--feature_preset",
                        dest="feature_preset",
                        action="store_true",
                        required=False,
                        help="Use the preset of features hardcoded.")

    return parser.parse_args()


def select_f_t(datadf, target_name, useClrt, ignore_col):
    try:
        targetv = datadf.pop(target_name).values
        labencode = LabelEncoder()

        targetv = labencode.fit_transform(targetv)  # Transform target column.

        datadf.drop(ignore_col,
                    axis="columns",
                    inplace=True)
    except KeyError:
        sys.exit("ERROR: No target found. Set the target column name!")
    if useClrt:
        featmtx = pd.DataFrame(clr(datadf))
    else:
        featmtx = datadf

    return targetv, featmtx

class RandomForestClassifierWithCoef(RFC):
    def fit(self, *args, **kwargs):
        super(RandomForestClassifierWithCoef, self).fit(*args, **kwargs)
        self.coef_ = self.feature_importances_


def select_features(X_train, y_train, cols, model):
    min_features_to_select = 1  # Minimum number of features to consider
    if model == "SVM":
        estimator = SVC(kernel="linear",
                        class_weight="balanced")
    elif model == "RF":
        estimator = RandomForestClassifierWithCoef(n_estimators=200,
                                                   bootstrap=True,
                                                   criterion="entropy",
                                                   max_depth=None,
                                                   max_features="auto",
                                                   min_samples_leaf=1,
                                                   min_samples_split=2,
                                                   class_weight="balanced")
    # if model == "KN":
    #     estimator = KNR(n_neighbors=5)
    # elif model == "LoR":
    #     estimator = LoR(solver="lbfgs",
    #                     max_iter=1000)
    # elif model == "SVM":
    #     estimator = SVC(kernel="linear")
    # elif model == "RF":
    #     estimator = RFC(criterion="entropy",
    #                     n_estimators=100)
    # elif model == "MLP":
    #     pass

    rfecv = RFECV(estimator=estimator,
                  step=1,
                  cv=StratifiedKFold(n_splits=9,
                                     shuffle=True),
                  scoring='f1',
                  min_features_to_select=min_features_to_select,
                  n_jobs=7)
    rfecv.fit(X_train, y_train)
    # Create RFE results.
    rfeDf = pd.DataFrame({"Features": cols[rfecv.support_],
                          "Ranking": rfecv.ranking_[rfecv.support_],
                          "Score": rfecv.grid_scores_[rfecv.support_],
                          "Importance": rfecv.estimator_.coef_[0]})

    # Plot number of features VS. cross-validation scores
    fig = plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Cross validation score (nb of correct classifications)")
    plt.plot(range(min_features_to_select,
                   len(rfecv.grid_scores_) + min_features_to_select),
             rfecv.grid_scores_)
    return(rfecv.get_support(1), rfeDf.sort_values(by=["Importance"]), fig)

def parameter_tuning(X_train, y_train):
    model_params = {"K-neighbors (KNR)": {"model": KNR(),
                                          "params": {'n_neighbors': [1, 2, 3, 4, 10]
                                                     }
                                    },
                    "Logistic_regression (LoR)": {"model": LoR(solver="lbfgs",
                                                               max_iter=1000,
                                                               class_weight="balanced"),
                                                  "params": {}},
                    # "Support_vector_machine (SVC)": {"model": SVC(probability=True,
                    #                                               class_weight="balanced"),
                    #                                  "params": {"C": [0.1, 1, 10, 100, 1000],
                    #                                             "gamma": [1, 0.1, 0.01, 0.001, 0.0001],
                    #                                             "kernel": ["rbf", "linear"]}
                    #                                  },
                    "Random_forest (RFC)": {"model": RFC(class_weight="balanced",
                                                         n_jobs=7,
                                                         verbose=0),
                                            "params": {"criterion": ["gini", "entropy"],
                                                       "n_estimators": [200],
                                                       "max_features": ["auto", 1, "sqrt", "log2"],
                                                       "max_depth": [None, 2, 5],
                                                       "min_samples_split": [2, 5],
                                                       "min_samples_leaf": [1, 2, 5],
                                                       "bootstrap": [True, False]}
                                            },
                    "Multi-layer_perceptron (MLP)": {"model": MLP(hidden_layer_sizes=(5,),
                                                                  random_state=42,
                                                                  max_iter=2000),
                                                     "params": {"solver": ["lbfgs", "sgd", "adam"],
                                                                "alpha": 10.0**-np.arange(1, 7)}
                                                     },
                    }

    scores = []
    best_score = 0
    print("=== HP tuning ===\n")
    for model_name, mod_par in model_params.items():
        print(model_name)

        grid = GCV(estimator=mod_par["model"],
                   param_grid=mod_par["params"],
                   cv=StratifiedKFold(n_splits=10,
                                      shuffle=True),
                                      n_jobs=7)
        grid.fit(X_train, y_train)
        scores.append({"model": model_name,
                       "best_score": grid.best_score_,
                       "best_params": grid.best_params_})

        if grid.best_score_ > best_score:
            best_model = grid  # Select model with the highest score.
            best_score = grid.best_score_

    print(pd.DataFrame(scores, columns=["model",
                                        "best_score",
                                        "best_params"]), "\n")
    print("Best model chosen:", best_model.best_estimator_)
    print("Best score:", best_model.best_score_)
    print("Best AUROC:", roc_auc_score(y_train,
                                       best_model.predict_proba(X_train)[:, 1]),
                                       "\n")
    sys.exit()
    return best_model


def train(X_train, y_train, mod_name, model, r_seed):
    # Selecting model with hyper parameter optimization followed.
    if mod_name == "PT":
        model_fit = model.fit(X_train, y_train)
    else:
        if mod_name == "KN":
            model = KNR(n_neighbors=5)
        elif mod_name == "LoR":
            model = LoR(solver="lbfgs",
                        max_iter=1000,
                        class_weight="balanced")
        elif mod_name == "RF":
            model = RFC(n_estimators=200,
                        bootstrap=True,
                        criterion="entropy",
                        max_depth=None,
                        max_features="auto",
                        min_samples_leaf=1,
                        min_samples_split=2,
                        class_weight="balanced")
        elif mod_name == "SVM":
            model = SVC(kernel="linear",
                        gamma=1,
                        C=100,
                        probability=True,
                        class_weight="balanced")
        elif mod_name == "MLP":
            model = MLP(hidden_layer_sizes=(5,),
                        random_state=42,
                        max_iter=2000,
                        solver="lbfgs",
                        alpha=0.0001)

        model_fit = model.fit(X_train, y_train)

    # ROC data for plotting.
    roc_Df = pd.DataFrame(list(zip(y_train,
                                   model.predict_proba(X_train)[:, 1])),
                          columns=["y_train",
                                   "probabilities"])
    roc_Df["r_seed"] = r_seed

    return model_fit, roc_Df


def predict(X_pre, y_pre, model_fit, r_seed):
    # Predicting using the fitted model.
    pred = model_fit.predict(X_pre)
    report = classification_report(y_pre,
                                   pred,
                                   target_names=["all",
                                                 "Cancer"],
                                   output_dict=True)
    auroc = roc_auc_score(y_pre,
                          model_fit.predict_proba(X_pre)[:, 1])

    # ROC data for plotting.
    fpr, tpr, thresholds = roc_curve(y_pre,
                                     model_fit.predict_proba(X_pre)[:, 1],
                                     pos_label=1)
    roc_Df = pd.DataFrame({"fpr": fpr,
                           "tpr": tpr,
                           "thresholds": thresholds})
    roc_Df["r_seed"] = r_seed

    tn, fp, fn, tp = confusion_matrix(y_pre,
                                      pred).ravel()
    rp = np.count_nonzero(y_pre == 1)  # Count real positives.
    rn = np.count_nonzero(y_pre == 0)  # Count real negatives.
    confmat_Df = pd.DataFrame({"true_positives": [tp],
                               "true_negatives": [tn],
                               "false_positives": [fp],
                               "false_negatives": [fn],
                               "r_seed": [r_seed]})

    inf = (tp/rp)-(fp/rn)  # informedness
    mk = (tp/(tp+fp))-(fn/(tn+fn)) # markedness
    ppv = tp/(tp+fp)
    npv = tn/(tn+fn)

    metricsDict = {"precision_1": report["Cancer"]["precision"],
                   "sensitivity": report["Cancer"]["recall"],
                   "specificity": report["all"]["recall"],
                   "accuracy": report["accuracy"],
                   "auroc": auroc,
                   "informedness": inf,
                   "markedness": mk,
                   "ppv": ppv,
                   "npv": npv,
                   "fdr": 1-ppv,
                   "for": 1-npv}

    preds_Df = pd.DataFrame({"Prediction": pred,
                             "Probability": model_fit.predict_proba(X_pre)[:, 1]}).set_index(X_pre.index)
    # preds_Df = pd.DataFrame()
    # print("=== Predicted confusion matrix === \n",
    #      confmat_Df)
    return metricsDict, roc_Df, confmat_Df, preds_Df


def main():
    args = argumentparsing()
    scaler = StandardScaler()

    if args.InputFile:
        datadf = pd.read_csv(args.InputFile,
                             sep=",")  # Read the dataset.

        datadf = datadf.sort_values(by=args.TargetColumn,
                                    ignore_index=True)
        samp_names = datadf[args.ignore_col]  # Separate sample names for later use.

        targetv, featmtx = select_f_t(datadf,
                                      args.TargetColumn,
                                      args.useClrt,
                                      args.ignore_col)
        # Split data.
        X_train, X_pre, y_train, y_pre = train_test_split(featmtx,
                                                          targetv,
                                                          stratify=targetv,
                                                          test_size=0.20,
                                                          random_state=42)
        # Scale the data using StandardScaler.
        scaler.fit(X_train)
        X_train = scaler.transform(X_train)
        X_pre = scaler.transform(X_pre)

        # Feature selection.
        if X_train.shape[1] > 3:
            selected_cols, rfeDf, fig = select_features(X_train,
                                                        y_train,
                                                        datadf.columns,
                                                        args.Models)
            print("\n=== Feature selection ===\n")
            print("Optimal number of features : %d" % len(selected_cols), "\n")

            # Save selected feature list and plot.
            rfeDf.to_csv("".join((args.OutputFolder,
                                  "_selected_features.csv")),
                                  index=False)
            fig.savefig("".join((args.OutputFolder, "_feature_selection.png")),
                        dpi=300)
        else:
            print("Feature count is < 3, feature selection is skipped!")
            selected_cols = list(range(0, len(featmtx.columns)))

        # parameter_tuning(X_train, y_train)
        sys.exit()
        if args.Models == "PT":
            # Model and param selection with selected features.
            best_model = parameter_tuning(X_train[:, selected_cols],
                                          y_train)
        else:
            best_model = "selected_model"

        roc_train_o = pd.DataFrame()
        roc_pre_o = pd.DataFrame()
        metrics_o = pd.DataFrame()
        confmat_Df_o = pd.DataFrame()
        inf_mark_o = pd.DataFrame()
        preds_Df_o = pd.DataFrame()
        for r_seed in [67, 98, 6, 14, 62, 12, 28, 13, 80, 79]:
            # Split data.
            X_train, X_pre, y_train, y_pre = train_test_split(featmtx,
                                                              targetv,
                                                              stratify=targetv,
                                                              test_size=0.20,
                                                              random_state=r_seed)
            # Model training.
            mod_fit, roc_train = train(X_train[X_train.columns[selected_cols]],
                                       y_train,
                                       args.Models,
                                       best_model,
                                       r_seed)
            roc_train_o = roc_train_o.append(roc_train,
                                             ignore_index=True)

            # Predictions.
            metricsDict, roc_pre, confmat_Df, preds_Df = predict(X_pre[X_pre.columns[selected_cols]],
                                                                 y_pre,
                                                                 mod_fit,
                                                                 r_seed)
            # Print the prediction score.
            print("Seed:", r_seed, "\n"
                  "Prediction score:", mod_fit.score(X_pre[X_pre.columns[selected_cols]],
                                                     y_pre))

            # Compute metrics.
            metricsDict_o = {"metrics": metricsDict,
                             "r_seed": r_seed}
            metrics_o = metrics_o.append(pd.DataFrame.from_dict({"value": metricsDict,
                                                                 "r_seed": r_seed}))
            metrics_o["metric"] = metrics_o.index

            roc_pre_o = roc_pre_o.append(roc_pre,
                                       ignore_index=True)

            confmat_Df_o = confmat_Df_o.append(confmat_Df,
                                               ignore_index=True)

            preds_Df[["WhichSample"]] = samp_names.loc[preds_Df.index]
            preds_Df[["r_seed"]] = r_seed
            preds_Df_o = preds_Df_o.append(preds_Df)
            preds_Df_o.sort_values(by="WhichSample",
                                   inplace=True)

        # Save metrics and ROC values.
        metrics_o.to_csv("".join((args.OutputFolder, "_metrics_predict.csv")),
                         index=False)
        roc_pre_o.to_csv("".join((args.OutputFolder, "_roc_predict.csv")),
                         index=False)  # Save roc file.
        roc_train_o.to_csv("".join((args.OutputFolder, "_roc_train.csv")),
                           index=False)
        confmat_Df_o.to_csv("".join((args.OutputFolder, "_confmat.csv")),
                            index=False)
        preds_Df_o.to_csv("".join((args.OutputFolder, "_preds.csv")),
                          columns=["WhichSample",
                                   "r_seed",
                                   "Prediction"],
                          index=False)

    elif (args.TrainingFile) and (args.PredictFile):
        cols = list(pd.read_csv(args.PredictFile, nrows=1))
        trainingDf = pd.read_csv(args.TrainingFile)  # Read the training data.
        trainingDf = trainingDf.sort_values(by=args.TargetColumn,
                                            ignore_index=True)

        predictDf = pd.read_csv(args.PredictFile)  # Read the prediction data.

        predictDf = predictDf.sort_values(by=args.TargetColumn,
                                          ignore_index=True)
        samp_names = predictDf[args.ignore_col]  # Separate sample names for later use.

        y_train, X_train = select_f_t(trainingDf,
                                      args.TargetColumn,
                                      args.useClrt,
                                      args.ignore_col)
        y_pre, X_pre = select_f_t(predictDf,
                                  args.TargetColumn,
                                  args.useClrt,
                                  args.ignore_col)

        # Feature selection.
        if args.feature_preset:
            selected_cols = [0,1,2,3,5,6,7,8,9,13,14,15,16,18,20,22,23,24,25,27,28,29,32,34,35,36,37,40,41,42,43,45,46,47,48,49,50,51,52,53,54,55,57,58,59,61,62,63]
        elif len(X_train.columns) > 3:
            selected_cols, rfeDf, fig = select_features(X_train,
                                                        y_train,
                                                        trainingDf.columns,
                                                        args.Models)
            print("\n=== Feature selection ===\n")
            print("Optimal number of features : %d" % len(selected_cols), "\n")
            # Save selected feature list and plot.
            rfeDf.to_csv("".join((args.OutputFolder,
                                  "_selected_features.csv")),
                                  index=False)
            fig.savefig("".join((args.OutputFolder, "_feature_selection.png")),
                        dpi=300)

        else:
            print("Feature count is < 3, feature selection is skipped!")
            selected_cols = list(range(0, len(X_train.columns)))

        # Scale the selected features using StandardScaler.
        scaler.fit(X_train[X_train.columns[selected_cols]])
        X_train = scaler.transform(X_train[X_train.columns[selected_cols]])
        X_pre = scaler.transform(X_pre[X_pre.columns[selected_cols]])
        X_pre = pd.DataFrame(X_pre)  # Convert to pdDf for indexing.

        parameter_tuning(X_train, y_train)

        roc_train_o = pd.DataFrame()
        roc_pre_o = pd.DataFrame()
        metrics_o = pd.DataFrame()
        confmat_Df_o = pd.DataFrame()
        inf_mark_o = pd.DataFrame()
        preds_Df_o = pd.DataFrame()
        mod_fit, roc_train = train(X_train,
                                   y_train,
                                   args.Models,
                                   "selected_model",
                                   42)  # Model training.
        metricsDict, roc_pre_o, confmat_Df_o, preds_Df = predict(X_pre,
                                                                 y_pre,
                                                                 mod_fit,
                                                                 42)  # Prediction.

        metrics_o = pd.DataFrame.from_dict({"value": metricsDict,
                                            "r_seed": 42})
        metrics_o["metric"] = metrics_o.index

        # Save metrics and ROC values.
        metrics_o.to_csv("".join((args.OutputFolder, "_metrics_predict.csv")),
                         index=False)
        roc_pre_o.to_csv("".join((args.OutputFolder, "_roc_predict.csv")),
                         index=False)  # Save roc file.
        roc_train_o.to_csv("".join((args.OutputFolder, "_roc_train.csv")),
                           index=False)
        confmat_Df_o.to_csv("".join((args.OutputFolder, "_confmat.csv")),
                            index=False)

        preds_Df[["WhichSample"]] = samp_names.loc[preds_Df.index]
        preds_Df[["r_seed"]] = 42
        preds_Df_o = preds_Df_o.append(preds_Df)
        preds_Df_o.sort_values(by="WhichSample",
                               inplace=True)
        preds_Df_o.to_csv("".join((args.OutputFolder, "_preds.csv")),
                          columns=["WhichSample",
                                   "r_seed",
                                   "Prediction",
                                   "Probability"],
                          index=False)

if __name__ == "__main__":
    main()
