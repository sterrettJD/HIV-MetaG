import os
import argparse
import pandas as pd


def get_args():
    """
    handles arg parsing for this script

    returns the parsed args
    """
    parser = argparse.ArgumentParser(
        prog="Aggregate Metaphlan bugs lists",
        description="Aggregates metaphlan bugs lists from humann outputs for multiple samples into a single tsv"
    )

    parser.add_argument("-i", "--indir",
                        required=True)
    parser.add_argument("-o", "--outfile",
                        required=True)

    parsed_args = parser.parse_args()
    return parsed_args


def get_filepaths(directory):
    # get list of only subdirectories
    subdirs = [d for d in os.listdir(directory) if os.path.isdir(d)]
    # assumes each of the subdirectories is named for a sample id,
    # and each of these contains a subdir sampleid.concat_humann_temp
    # and each of the bugs list is named sampleid.concat_metaphlan_bugs_list.tsv
    filepaths = [os.path.join(sampid,
                              f"{sampid}.concat_humann_temp",
                              f"{sampid}.concat_metaphlan_bugs_list.tsv") for sampid in subdirs]
    return filepaths, subdirs


def concat_files(filepaths_list, sampids_list):
    for i, f in enumerate(filepaths_list):
        sample_id = sampids_list[i]

        df = pd.read_csv(f, sep="\t", header=4)
        # set index to be the clade name so we can merge the dataframes
        df.index = df["#clade_name"]
        df.rename({"relative_abundance": sample_id})
        # create our full df if it doesn't exist yet
        if i == 0:
            full = df.copy()
        # add the new sample to the full df
        full = pd.concat(full, df[sample_id], axis=1)

    return full

if __name__ == "__main__":
    args = get_args()

    filepaths_list, sampids_list = get_filepaths(args.indir)
    full = concat_files(filepaths_list, sampids_list)

    full.to_csv(f"{args.outfile}", sep="\t")
