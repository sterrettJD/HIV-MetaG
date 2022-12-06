import pandas as pd
import numpy as np
import argparse
import os


def get_args():
    """
    handles arg parsing for this script

    returns the args
    """
    parser = argparse.ArgumentParser(
        prog="Demux to metadata",
        description="Converts Anschutz Genomics Core demux files to YMP metadata"
    )

    parser.add_argument("-i", "--infile",
                        required=True)
    parser.add_argument("-o", "--outfile",
                        required=True)
    parser.add_argument("-d", "--directory",
                        required=True)

    parsed_args = parser.parse_args()
    return parsed_args


def sample_id_to_files(dir, sampleid, read_sig):
    """
    :param dir: str: directory to search
    :param sampleid: str: sampleid to look for
    :param read_sig: str: signifier to also look for (e.g., "R1" for fwd read)
    :return: str: filename
    """
    for file in os.listdir(dir):
        if file.startswith(sampleid) and read_sig in file:
            return file

    # if no file returned
    return np.nan


if __name__ == "__main__":

    args = get_args()

    df = pd.read_csv(args.infile,
                     skiprows=5).drop(0, axis=0)

    df["ForwardReads"] = df["Sample"].apply(
        lambda x: sample_id_to_files(args.directory, x, "R1")
    )
    df["ReverseReads"] = df["Sample"].apply(
        lambda x: sample_id_to_files(args.directory, x, "R2")
    )

    order = ["Sample", "ForwardReads", "ReverseReads"]
    order = order + [col for col in df.columns if col not in order]

    new_df = df[order]

    print([col for col in new_df.columns])
    print(new_df.head())

    new_df.to_csv(args.outfile, index=False)
