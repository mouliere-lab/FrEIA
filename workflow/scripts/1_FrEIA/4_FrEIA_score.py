#!/usr/bin/env python3
# Author: Norbert Moldován
# Mouliere Lab
# Amsterdam UMC
# Vrije Universiteit Amsterdam
import sys
from multiprocessing import Pool
from functools import partial
import argparse
from FrEIA_tools import CastDataTypes
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.spatial import distance
from combat.pycombat import pycombat
import umap
import matplotlib
matplotlib.use('Agg')  # Matplotlib can't use interactive backend on HPC.
import matplotlib.pyplot as plt
import seaborn as sns

# This scritp will calculate the FrEIA score of samples
# by using trinucleotide fragment end sequence proportions and diversity scores.
# Batch correction can also be applied on the fragment end sequence proportions
# using pycombat.


def argumentparsing():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_path",
                        dest="input",
                        type=str,
                        required=True,
                        help="The path to the FrEIA directory "
                        "with trinucleotide FrES proportions and diversity.")
    parser.add_argument("-o", "--output",
                        dest="outpath",
                        type=str,
                        required=True,
                        help="The path to the output file.")
    parser.add_argument("-p", "--median_panel",
                        dest="median_panel",
                        type=str,
                        required=False,
                        help="The path to the file containing the median panel."
                             "This can be generated using the Create_panels.py.")
    parser.add_argument("-m", "--metadata_path",
                        dest="metadata",
                        type=str,
                        required=True,
                        help="The path to the metadata.csv file.")
    parser.add_argument("--control_name",
                        dest="ctr_name",
                        type=str,
                        required=True,
                        help="How controls are named in the metadata file.")
    parser.add_argument("--case_name",
                        dest="cas_name",
                        type=str,
                        required=True,
                        help="How cases are named in the metadata file.")
    parser.add_argument("--tnc_control",
                        dest="tnc_ctr_in",
                        type=str,
                        required=False,
                        help="A coma separated list of trinucleotides to use when computing the FrEIA score."
                             "If not given it will be calculated from the data as the trinucleotides,"
                             "increased in controls.")
    parser.add_argument("--tnc_case",
                        dest="tnc_cas_in",
                        type=str,
                        required=False,
                        help="A coma separated list of trinucleotides to use when computing the FrEIA score."
                             "If not given it will be calculated from the data as the trinucleotides,"
                             "increased in cases.")
    parser.add_argument("--gini",
                        dest="use_gini",
                        action="store_true",
                        required=False,
                        help="If set the Gini indexes are used in "
                        "the FrEIA score calculation.")
    parser.add_argument("--mds",
                        dest="use_mds",
                        action="store_true",
                        required=False,
                        help="If set the normalized Shannon entropy is "
                        "used in the FrEIA score calculation.")
    parser.add_argument("-b",
                        dest="batches",
                        type=str,
                        required=False,
                        help="If given batch correction will be performed using "
                        "the provided variable. The variable has to be a column "
                        "name in the metadata file.")
    parser.add_argument("-t", "--threads",
                        dest="threads",
                        type=int,
                        required=True,
                        help="Number of threads.")
    return parser.parse_args()


def mannwhitneyu_test_worker(dat_Df, args, b):
    dat_Df = dat_Df[dat_Df.base == b]

    if (len(dat_Df[dat_Df.phenotype == args.cas_name].Proportion) > 0 and
       len(dat_Df[dat_Df.phenotype == args.ctr_name].Proportion) > 0):
        U1, p = mannwhitneyu(dat_Df[dat_Df.phenotype == args.cas_name].Proportion.tolist(),
                             dat_Df[dat_Df.phenotype == args.ctr_name].Proportion.tolist())
    else:
        p = 1

    out_Df = pd.DataFrame({"base": [b],
                           "p-val": [p]})
    return out_Df


def mannwhitneyu_test(dat_Df, meta_df, args):

    dat_Df = dat_Df.merge(meta_df[["sample_name",
                                   "phenotype"]])
    # Compare the distribution of phenotypes per base and
    # compute the p-values per base.
    mwu_Df = pd.DataFrame()
    base_L = dat_Df.base.unique().tolist()

    pool = Pool(processes=args.threads)

    mwu_Df = mwu_Df.append(pool.map(partial(mannwhitneyu_test_worker,
                                            dat_Df,
                                            args),
                                    base_L))
    pool.close()
    pool.join()

    mwu_Df.reset_index(drop=True,
                       inplace=True)
    return mwu_Df


def logfc(dat_Df, meta_df, args):
    dat_Df = dat_Df.merge(meta_df[["sample_name",
                                   "phenotype"]])
    # 1. Calculate the median proportion per base.
    median_Df = dat_Df.groupby(["base", "phenotype"],
                               observed=True).Proportion.median().reset_index(level="phenotype")

    # 2. Compute log10FC per base.
    control_median_L = median_Df[median_Df.phenotype == args.ctr_name].Proportion
    case_median_L = median_Df[median_Df.phenotype == args.cas_name].Proportion

    out_Df = np.log10(case_median_L/control_median_L).reset_index()
    out_Df.rename(columns={"Proportion": "Log10FC"},
                  inplace=True)

    return out_Df, control_median_L, case_median_L


def get_signif_bases(dat_Df):
    # 1. Compute the 5% and 95% quantiles of log10FC.
    Q25 = dat_Df.Log10FC.quantile(q=0.25)
    Q75 = dat_Df.Log10FC.quantile(q=0.75)
    # 2. Subset base with a p-val < 0.05 and outside of the quantile range.
    control_base_L = dat_Df[(dat_Df["p-val"] < 0.01) &
                           (dat_Df["Log10FC"] < Q25)].base.tolist()
    case_base_L = dat_Df[(dat_Df["p-val"] < 0.01) &
                          (dat_Df["Log10FC"] > Q75)].base.tolist()
    #print(case_base_L)
    return control_base_L, case_base_L


def umap_projection(dat_Df, batch_L, val_name, args):
    print("\n === UMAP projection of controls === \n")
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(dat_Df.T)
    plt_Df = pd.DataFrame()
    plt_Df["PC0"] = embedding[:, 0]
    plt_Df["PC1"] = embedding[:, 1]
    plt_Df["category"] = batch_L
    # Plot UMAP projection.
    sns.scatterplot(x="PC0",
                    y="PC1",
                    data=plt_Df,
                    hue="category")
    plt.gca().set_aspect('equal',
                         'datalim')
    plt.title(" ".join(('UMAP projection of',
                        val_name,
                        'control samples.')),
              fontsize=11)
    plt.savefig("".join((args.outpath,
                         "_controls_UMAP_",
                         val_name,
                         ".png")))
    plt.close()


def batch_correction(dat_Df, batch_L, control_L, val_name, args):
    # Reshape the Df and remove NaN-s.
    if val_name == "Proportion":
        dat_Df.drop(["group"],
                    axis=1,
                    inplace=True)
    else:
        dat_Df = dat_Df[["sample_name",
                         val_name]]

    dat_Df = dat_Df.transpose().dropna()
    dat_Df.columns = dat_Df.loc["sample_name"]
    dat_Df.drop(["sample_name"],
                inplace=True)

    if args.batches != 'None':  # Perform batch correction only if requested!
        # UMAP projection of controls.
        umap_projection(dat_Df[control_L],
                        batch_L[batch_L.phenotype == args.ctr_name][args.batches].reset_index(drop=True),
                        "_".join(("raw",
                                  val_name)),
                        args)

        # Perform batch correction.
        out_Df = pycombat(dat_Df.astype("float"), batch_L[args.batches])

        try:
            # UMAP projection.
            umap_projection(out_Df[control_L],
                            batch_L[batch_L.phenotype == args.ctr_name][args.batches].reset_index(drop=True),
                            "_".join(("corrected",
                                      val_name)),
                            args)
        except ValueError:
            print("Batch correction resulted in NaNs. UMAP projection skipped!")

    else:  # If batch correction is not requested return the transposed dataframe.
        out_Df = dat_Df

    out_Df["base"] = out_Df.index

    # Reshape the Df to long.
    out_Df = out_Df.melt(id_vars=["base"],
                         var_name="sample_name",
                         value_name=val_name)
    return out_Df


def euclidean_worker(samp_df, control_df, case_df, samp):
    # Select sample for worker.
    samp_df = samp_df[samp_df.sample_name == samp]

    # Subset selected bases to those that are present in all samples.
    control_df = samp_df.merge(control_df,
                               on="base",
                               suffixes=["_sample",
                                         "_control"])
    case_df = samp_df.merge(case_df,
                              on="base",
                              suffixes=["_sample",
                                        "_case"])
    control_dist = distance.euclidean(control_df.Proportion_sample,
                                      control_df.Proportion_control)
    case_dist = distance.euclidean(case_df.Proportion_sample,
                                     case_df.Proportion_case)
    FrEIA_score = np.sqrt((control_dist**2)/(case_dist**2))

    out_Df = pd.DataFrame({"sample_name": [samp],
                           "control_ED": [control_dist],
                           "case_ED": [case_dist],
                           "FrEIA_score": [FrEIA_score]})

    return out_Df


def euclidean_FrEIA_score(prop_df, div_df, control_V, case_V, args):
    control_df = control_V.reset_index()
    case_df = case_V.reset_index()

    if args.use_gini:
        # Calculate the median diversity per bases for phenotype.
        median_div_Df = div_df.groupby(["base", "phenotype"],
                                       observed=True).Gini.median().reset_index(level=["base",
                                                                                       "phenotype"])
        # Concat diversity score to proportions.
        div_df.rename(columns={"Gini": "Proportion"},
                      inplace=True)
        median_div_Df.rename(columns={"Gini": "Proportion"},
                      inplace=True)
        samp_df = pd.concat([prop_df[["sample_name", "base", "Proportion"]],
                             div_df[["sample_name", "base", "Proportion"]]])
        control_df = pd.concat([control_df[["base", "Proportion"]],
                                median_div_Df[median_div_Df.phenotype == args.ctr_name][["base", "Proportion"]]])
    elif args.use_mds:
        # Calculate the median diversity per bases for phenotype.
        median_div_Df = div_df.groupby(["base", "phenotype"],
                                       observed=True).MDS.median().reset_index(level=["base",
                                                                                      "phenotype"])
        # Concat diversity score to proportions.
        div_df.rename(columns={"MDS": "Proportion"},
                      inplace=True)
        median_div_Df.rename(columns={"MDS": "Proportion"},
                             inplace=True)
        samp_df = pd.concat([prop_df[["sample_name", "base", "Proportion"]],
                             div_df[["sample_name", "base", "Proportion"]]])
        case_df = pd.concat([case_df[["base", "Proportion"]],
                               median_div_Df[median_div_Df.phenotype == args.cas_name][["base", "Proportion"]]])
    else:
        samp_df = prop_df[["sample_name", "base", "Proportion"]]

    # Compute the distance of each sample from the control and
    # case median vectors.
    ed_Df = pd.DataFrame()
    pool = Pool(processes=args.threads)

    samp_L = prop_df.sample_name.unique().tolist()
    ed_Df = ed_Df.append(pool.map(partial(euclidean_worker,
                                          samp_df,
                                          control_df,
                                          case_df),
                                  samp_L))
    pool.close()
    pool.join()

    return ed_Df


def main():
    args = argumentparsing()  # Parse arguments.

    prop_Df = pd.read_csv("".join((args.input, "/", "Dat_GT__sample.csv")),
                          sep=",")  # Read proportions.
    div_Df = pd.read_csv("".join((args.input, "/", "Dat_GT__MDS_sample.csv")),
                          sep=",")  # Read diversity.
    meta_Df = pd.read_csv(args.metadata,
                          sep=",")  # Read metadata.
    if args.median_panel:
        med_pan_Df = pd.read_csv(args.median_panel,
                                 sep=",")  # Read panel of medians.

    prop_Df.rename(columns={"WhichSample": "sample_name",
                            "WhichGroup": "group"},
                   inplace=True)
    prop_Df = CastDataTypes(prop_Df)  # Cast to cathegorycal.
    div_Df = CastDataTypes(div_Df)  # Cast to cathegorycal.

    control_L = meta_Df[meta_Df.phenotype == args.ctr_name].sample_name.tolist()

    if args.batches != 'None':  # Perform batch correction.
        batch_L = meta_Df[["phenotype", args.batches]]
    else:
        batch_L = ["None"]

    prop_Df = batch_correction(prop_Df, batch_L, control_L, "Proportion", args)
    prop_Df.Proportion = pd.to_numeric(prop_Df.Proportion)

    if args.use_mds:
        div_Df = batch_correction(div_Df, batch_L, control_L, "MDS", args)
        div_Df.MDS = pd.to_numeric(div_Df.MDS)
    elif args.use_gini:
        div_Df = batch_correction(div_Df, batch_L, control_L, "Gini", args)
        div_Df.Gini = pd.to_numeric(div_Df.Gini)

    if args.batches != 'None':
        # Save the corrected data.
        prop_Df.to_csv("".join((args.outpath,
                                "_corrected_tnc.csv")),
                       index=False)
        div_Df.to_csv("".join((args.outpath,
                               "_corrected_div.csv")),
                      index=False)

    mwu_Df = mannwhitneyu_test(prop_Df,
                               meta_Df,
                               args)  # Compare distributions per base.
    logfc_Df, control_median_L, case_median_L = logfc(prop_Df,
                                                      meta_Df,
                                                      args)  # Compute log10FC per base.
    base_out_Df = logfc_Df.merge(mwu_Df)
    base_out_Df.to_csv("".join((args.outpath,
                                "_logFC_p.csv")),
                       index=False)

    if args.median_panel:  # If panel of medians is provided use that instead of the generated one.
        control_median_L = med_pan_Df[["base", "ctr"]].rename(columns={"ctr": "Proportion"})
        control_median_L.set_index("base",
                                   inplace=True)
        control_median_L = control_median_L.squeeze()

        case_median_L = med_pan_Df[["base", "case"]].rename(columns={"case": "Proportion"})
        case_median_L.set_index("base",
                                inplace=True)
        case_median_L = case_median_L.squeeze()

    if args.tnc_ctr_in and args.tnc_cas_in:
        control_base_L = args.tnc_ctr_in.split(",")
        case_base_L = args.tnc_cas_in.split(",")
    else:
        # Compute significant bases.
        control_base_L, case_base_L = get_signif_bases(base_out_Df)

    # Subset prop_Df to contain only the significant bases.
    prop_Df = prop_Df[prop_Df.base.isin(control_base_L + case_base_L)].reset_index(drop=True)

    # Add the phenotype column to diversity Df.
    div_Df = div_Df.merge(meta_Df[["sample_name",
                                   "phenotype"]])
        # Compute the FrEIA score.
    FrEIA_score_df = euclidean_FrEIA_score(prop_Df,
                                           div_Df,
                                           control_median_L.loc[control_base_L],
                                           case_median_L.loc[case_base_L],
                                           args)
    FrEIA_score_df.to_csv("".join((args.outpath,
                                   "_FrEIA_score.csv")),
                          index=False)


if __name__ == "__main__":
    main()
