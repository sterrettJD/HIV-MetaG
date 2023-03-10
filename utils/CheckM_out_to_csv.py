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
    parser.add_argument("-i", "--in",
                        required=True)
    parser.add_argument("-o", "--out",
                        required=True)

    parsed_args = parser.parse_args()
    return parsed_args


def main():
    args = get_args()
    # weird delimiter of at least 2 spaces seems to be working. It's not a tab, and some entries contain spaces
    df = pd.read_csv(args.in, skip_footer=2, delimiter=r"\s\s+")

    df.drop(0, inplace=True)
    df.to_csv(args.out)

if __name__=="__main__":
    main()
