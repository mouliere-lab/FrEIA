#!/usr/bin/env python3
# Norbert Moldován

import os
import sys
import subprocess
from multiprocessing import Pool
from functools import partial
import argparse
from itertools import product
import pandas as pd
from FrEIA_tools import CastDataTypes


def parsingarguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fregmentEndFile",
                        dest="fregmentEndFile",
                        type=str,
                        required=True,
                        help="The path to the sample file "
                        "produced in extraction step."
                        "This has to be a 'Fragment ends' file!")
    parser.add_argument("-bam", "--sampleBam",
                        dest="SampleBam",
                        type=str,
                        required=True,
                        help="The path to the bam file.")
    parser.add_argument("-o", "--output",
                        dest="outpath",
                        type=str,
                        required=True,
                        help="The path to the outputput folder.")
    parser.add_argument("-st", "--subsetByFeature",
                        dest="subsetByFeature",
                        type=str,
                        help="The path to a list of features present in the "
                        "'Fragment ends' file to subset the bam. "
                        "If this is given subsetting by trinucleotides "
                        "will be turned off!")
    parser.add_argument("-t", "--threads",
                        dest="threads",
                        type=int,
                        required=True,
                        help="Number of threads.")
    return parser.parse_args()


def readdata(filepath):
    try:
        dataframe = pd.read_parquet(filepath)
    # dataframe = dataframe.sample(n=100000)  # subsampling to evade killing
    except FileNotFoundError:
        sys.exit("No datafiles found! Please supply the appropriate files! \n")

    try:
        dataframe = CastDataTypes(dataframe)
    except KeyError:
        sys.exit("Error during data type casting! \n")
    return dataframe


def subsetreadlist(dataframe, filename, outpath, feature):
    featuresubset = dataframe[dataframe["leftmost_tnc_seq"] == feature]

    featuresubset["qname"].to_csv("".join((outpath,
                                           "/subset_read_list/feature/",
                                           filename, "_",
                                           feature, ".txt")),
                                  index=False,
                                  header=False)  # Subset by feature.

    feature_readc = len(featuresubset)
    readcsubset = dataframe.sample(n=feature_readc,
                                   replace=True,
                                   random_state=feature_readc
                                   ).reset_index(drop=True)
    readcsubset["qname"].to_csv("".join((outpath,
                                         "/subset_read_list/readcount/",
                                         filename, "_",
                                         feature, ".txt")),
                                index=False,
                                header=False)  # Subset by readcount.


def subsetbam(bampath, outpath, filename, whichsubset_feature):
    whichsubset, feature = whichsubset_feature
    PicardFilterReads = "picard FilterSamReads \
                        I={} \
                        O={} \
                        READ_LIST_FILE={} \
                        FILTER=includeReadList".format(
        bampath,
        "".join((outpath, "/subset_bams/",
                 whichsubset, "/",
                 filename, "_",
                 feature, ".bam")),
        "".join((outpath, "/subset_read_list/",
                 whichsubset, "/",
                 filename, "_",
                 feature, ".txt")))
    subprocess.run(PicardFilterReads, shell=True)


def createfolders(outpath):
    if not os.path.exists("".join((outpath, "/subset_read_list/feature/"))):
        os.makedirs("".join((outpath, "/subset_read_list/")))
        os.makedirs("".join((outpath, "/subset_read_list/feature/")))
        os.makedirs("".join((outpath, "/subset_read_list/readcount/")))
    if not os.path.exists("".join((outpath, "/subset_bams/feature/"))):
        os.makedirs("".join((outpath, "/subset_bams/")))
        os.makedirs("".join((outpath, "/subset_bams/feature/")))
        os.makedirs("".join((outpath, "/subset_bams/readcount/")))


def main():
    args = parsingarguments()
    dataframe = readdata(args.fregmentEndFile)

    createfolders(args.outpath)

    filename = args.fregmentEndFile.split("/")[-1].strip(".pq")
    featureset = set(dataframe[~dataframe.leftmost_read_seq.str.contains("N")]["leftmost_tnc_seq"])

    whichsubset_feature = list(product(["feature", "readcount"], featureset))

    pool = Pool(processes=args.threads)

    pool.map(partial(subsetreadlist,
                     dataframe,
                     filename,
                     args.outpath), featureset)

    pool.map(partial(subsetbam,
                     args.SampleBam,
                     args.outpath,
                     filename), whichsubset_feature)

    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
