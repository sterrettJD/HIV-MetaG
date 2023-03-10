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


def main():
    args = get_args()
    # weird delimiter of at least 2 spaces seems to be working. It's not a tab, and some entries contain spaces
    df = pd.read_csv(args.infile, skipfooter=2, delimiter=r"\s\s+", engine="python")

    df.drop(0, inplace=True)
    df.to_csv(args.outfile, index=False)

if __name__=="__main__":
    main()
