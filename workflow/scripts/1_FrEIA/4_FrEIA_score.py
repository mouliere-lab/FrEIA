#!/usr/bin/env python3
# Author: Norbert Moldován
import os
import sys
from multiprocessing import Pool
from functools import partial
import argparse
from FrEIA_tools import (CastDataTypes)
import pandas as pd


def argumentparsing():
    parser = argparse.ArgumentParser()
    parser.add_argument("-ip", "--input_proportions",
                        dest="in_prop",
                        type=str,
                        required=True,
                        help="The path "
                        "to the files containing FrES proportions.")
    parser.add_argument("-id", "--input_diversity",
                        dest="in_div",
                        type=str,
                        required=True,
                        help="The path "
                        "to the files containing FrES diversities.")
    parser.add_argument("-st", "--sample_table",
                        dest="samp_tb",
                        type=str,
                        required=True,
                        help="The path "
                        "to the sample table.")
    parser.add_argument("-o", "--output",
                        dest="outpath",
                        type=str,
                        required=True,
                        help="The path "
                        "to the output file.")

    return parser.parse_args()

def select_trinucs():  # Select trinucleotides passing the threshold.
    trinc_high = ["AAG",
                  "ACG",
                  "AGA",
                  "AGG",
                  "ATA",
                  "ATG",
                  "GAG",
                  "GCG",
                  "GTA",
                  "GTG",
                  "TCG",
                  "TTA",
                  "TTC",
                  "TTG"]
    trinc_low = ["ACA",
                 "ACC",
                 "ACT",
                 "CAC",
                 "CAT",
                 "CCA",
                 "CCC",
                 "CCT",
                 "CGA",
                 "CGC",
                 "CGG",
                 "CGT",
                 "CTT",
                 "GCC",
                 "GGC"]
    return trinc_high, trinc_low

def FrEIA_score(stdf, propdf, divdf, trinc):  # Calculates the FrEIA score.
    trinc_high, trinc_low = trinc
    FrEIADf = pd.DataFrame()
    for sample in stdf.sample_name:
        dd = 0
        dr = 0
        for tnc in trinc_high:
            dd = dd + propdf.loc[propdf.sample_name == sample,
                                 tnc].values  # Sum of proportions of trinc higher in cancer.
        for tnc in trinc_low:
            dr = dr + propdf.loc[propdf.sample_name == sample,
                                 tnc].values  # Sum of proportions of trinc lower in cancer.
        S = divdf.loc[divdf.sample_name == sample,
                      "MDS"].values  # Shannon entropy of sample.
        FrEIA_score = (dd/dr)*S  # FrEIA score.
        # print(FrEIA_score)
        FDf = pd.DataFrame(data={"sample_name": sample,
                                 "FrEIA_score": FrEIA_score})
        FrEIADf = FrEIADf.append(FDf)
        FrEIADf = FrEIADf.reset_index(drop=True)
    return FrEIADf

def main():
    args = argumentparsing()  # Parse arguments.

    stdf = pd.read_csv(args.samp_tb, sep=" ")
    propdf = pd.read_csv(args.in_prop)
    propdf = propdf.fillna(propdf.groupby(["WhichGroup"]).median())
    propdf = propdf.rename(columns={"WhichSample": "sample_name",
                                    "WhichGroup": "group"})
    propdf = pd.merge(propdf, stdf)

    divdf = pd.read_csv(args.in_div)
    divdf = divdf.fillna(divdf.groupby(["group"]).median())
    divdf = pd.merge(divdf, stdf)

    FrEIADf = FrEIA_score(stdf, propdf, divdf, select_trinucs())
    FrEIADf.to_csv(args.outpath,
                   index=False)

if __name__ == "__main__":
    main()
