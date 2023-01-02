import argparse
from aggregate_metaphlan_bugslists import get_filepaths


def get_args():
    """
    handles arg parsing for this script

    returns the parsed args
    """
    parser = argparse.ArgumentParser(
        prog="Convert MetaPhlan v4 reports to v3 format",
        description="Convert MetaPhlan v4 reports to v3 format"
    )

    parser.add_argument("-i", "--indir",
                        required=True)

    parsed_args = parser.parse_args()
    return parsed_args


def fix_files(filepaths):
    for file in filepaths:
        with open(file, "r") as f:
            contents = f.readlines()

        # replace the first line to pass the Pavian checks for Metaphlan v3
        contents[0] = "mpa_v3"

        # create new filepath
        # should take "dir/subdir/bugs_list.tsv" -> "dir/subdir/bugs_list_v3.tsv"
        split_by_dot = file.split(".")
        split_by_dot[-2] = split_by_dot[-2] + "_v3"
        v3_path = ".".join(split_by_dot)

        # write it to a new file
        with open(v3_path, "w") as outfile:
            outfile.writelines(contents)


if __name__ == "__main__":
    args = get_args()
    filepaths = get_filepaths(args.indir)
    fix_files(filepaths)
