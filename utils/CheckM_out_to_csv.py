import argparse
import pandas as pd

def get_args():
    """
    handles arg parsing for this script

    returns the parsed args
    """
    parser = argparse.ArgumentParser(
        prog="CheckM Output Converter",
        description="Parses CheckM lineage_wf output and converts it to a usable csv"
    )
    parser.add_argument("-i", "--infile",
                        required=True)
    parser.add_argument("-o", "--outfile",
                        required=True)

    parsed_args = parser.parse_args()
    return parsed_args

def replace_bin_id_with_filepath(df):
    new_df = df.copy()
    new_df["Bin Id"] = new_df["Bin Id"].apply(lambda x: "../hiv.t32.n40.metaspades.metabat2/bins_to_derep/" + x + ".fa")
    return new_df

def rename_columns(df):
    return df.rename(columns={"Bin Id": "genome",
                              "Completeness": "completeness",
                              "Contamination": "contamination",
                              "Strain heterogeneity": "heterogeneity"})

def main():
    args = get_args()
    # weird delimiter of at least 2 spaces seems to be working. It's not a tab, and some entries contain spaces
    df = pd.read_csv(args.infile, skipfooter=2, delimiter=r"\s\s+", engine="python")
    df.drop(0, inplace=True)
    df = replace_bin_id_with_filepath(df)
    df = rename_columns(df)
    df.to_csv(args.outfile, index=False)

if __name__=="__main__":
    main()
