#!/usr/bin/env python3
# Norbert Moldovan

import os
import sys
import subprocess
from multiprocessing import Pool
from functools import partial
import argparse
from itertools import product
import pandas as pd
from FrEIA_tools import CastDataTypes


# Parsing the arguments + some nice help text.
def parsingarguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-iFE", "--inputFE",
                        dest="inPathFE",
                        type=str,
                        required=True,
                        help="The "
                        "path to the sample file to be diluted."
                        "This has to be a 'Fragment ends' file!")
    parser.add_argument("-sFE", "--solventFE",
                        dest="solventPathFE",
                        type=str,
                        required=True,
                        help="The "
                        "path to the file to dilute with."
                        "This has to be a 'Fragment ends' file!")
    parser.add_argument("-iB", "--inputB",
                        dest="inPathB",
                        type=str,
                        required=False,
                        help="The "
                        "path to the sample file to be diluted."
                        "This has to be a bam file!")
    parser.add_argument("-sB", "--solventB",
                        dest="solventPathB",
                        type=str,
                        required=False,
                        help="The "
                        "path to the file to dilute with."
                        "This has to be a bam file!")
    parser.add_argument("-o", "--output",
                        dest="outPath",
                        type=str,
                        required=True,
                        help="Output file prefix with path.")
    parser.add_argument("-ic", "--inputConcentration",
                        dest="inCon",
                        type=float,
                        required=True,
                        help="The concentration of the input as fraction.")
    parser.add_argument("-r", "--replicate",
                        dest="replicate",
                        type=int,
                        required=False,
                        help="The number of replication. Default is 1.")
    parser.add_argument("-t", "--threads",
                        dest="threads",
                        type=int,
                        required=True,
                        help="Number of threads.")
    return parser.parse_args()


def readdata(filepath):
    try:
        dataframe = pd.read_parquet(filepath)
        # dataframe = dataframe.sample(n=1000000)  # subsampling to evade killing
    except FileNotFoundError:
        sys.exit("No datafiles found! Please supply the appropriate files! \n")

    try:
        dataframe = CastDataTypes(dataframe)
    except KeyError:
        sys.exit("Error during data type casting! \n")
    return dataframe


def dilutepq(solute, solvent, incon, outpath, finalcon_r_src):
    finalcon, r, src = finalcon_r_src

    final_rc = len(solvent)
    dilutionfactor = incon / finalcon
    solute_rc = int(final_rc // dilutionfactor)
    # print(final_rc, dilutionfactor, solute_rc)
    solutedf = solute.sample(n=solute_rc,
                             replace=True,
                             random_state=solute_rc*r)
    solventdf = solvent.sample(n=final_rc - solute_rc,
                               replace=True,
                               random_state=solute_rc*r)

    dilution = solventdf.append(solutedf)

    # Subsample dataset to given size.
    dilution = dilution.sample(n=src,
                               random_state=solute_rc*r).reset_index(drop=True)

    dilution.to_parquet("".join((outpath, "/diluted_fragment_ends/"
                                 "TF_", str(finalcon),
                                 "_RC_", str(src),
                                 "_iter_", str(r + 1),
                                 ".pq")),
                        compression="gzip",
                        index=False)  # Save diluted 'Fragment ends' file.

    dilution["qname"].to_csv("".join((outpath, "/diluted_read_lists/",
                                      "TF_", str(finalcon),
                                      "_RC_", str(src),
                                      "_iter_", str(r + 1),
                                      ".txt")),
                             index=False)  # Save readname file.

    return dilution


def mergebam(solute, solvent, outpath, threads):

    SamtoolsMerge = "samtools merge -@ {} {} {} {}".format(
                    threads,
                    "".join((outpath, "/solute_solvent_merged.bam")),
                    solute,
                    solvent)
    subprocess.run(SamtoolsMerge, shell=True)


def dilutebam(outpath, finalcon_r_src):
    finalcon, r, src = finalcon_r_src

    PicardFilterReads = "picard FilterSamReads \
                        I={} \
                        O={} \
                        READ_LIST_FILE={} \
                        FILTER=includeReadList".format(
        "".join((outpath, "/solute_solvent_merged.bam")),
        "".join((outpath, "/diluted_bams/"
                 "TF_", str(finalcon),
                 "_RC_", str(src),
                 "_iter_", str(r + 1),
                 ".bam")),
        "".join((outpath, "/diluted_read_lists/",
                 "TF_", str(finalcon),
                 "_RC_", str(src),
                 "_iter_", str(r + 1),
                 ".txt")))
    subprocess.run(PicardFilterReads, shell=True)


def createfolders(outpath):
    if not os.path.exists("".join((outpath, "/diluted_fragment_ends/"))):
        os.makedirs("".join((outpath, "/diluted_fragment_ends/")))
    if not os.path.exists("".join((outpath, "/diluted_read_lists/"))):
        os.makedirs("".join((outpath, "/diluted_read_lists/")))
    if not os.path.exists("".join((outpath, "/diluted_bams/"))):
        os.makedirs("".join((outpath, "/diluted_bams/")))


def main():
    args = parsingarguments()
    solute = readdata(args.inPathFE)
    solvent = readdata(args.solventPathFE)

    createfolders(args.outPath)

    solutionrc = [1000000, 500000, 100000,
                  50000, 10000, 1000]  # Used to downsample the final dataset.

    finalcon = [0.001, 0.0025, 0.005, 0.0075,
                0.01, 0.025, 0.05, 0.075, 0.1]  # The acheavable concentration.

    if not args.replicate:
        replist = [1]
    else:
        replist = list(range(args.replicate))

    finalcon_r_src = list(product(finalcon, replist, solutionrc))

    mergebam(args.inPathB,
             args.solventPathB,
             args.outPath,
             args.threads)

    pool = Pool(processes=args.threads)
    pool.map(partial(dilutepq,
                     solute,
                     solvent,
                     args.inCon,
                     args.outPath), finalcon_r_src)

    pool.map(partial(dilutebam,
                     args.outPath), finalcon_r_src)

    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
