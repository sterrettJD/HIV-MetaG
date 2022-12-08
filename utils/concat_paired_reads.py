import subprocess
import os
import argparse
import numpy as np

def get_args():
    """
    handles arg parsing for this script

    returns the parsed args
    """
    parser = argparse.ArgumentParser(
        prog="Concat paired reads",
        description="Concatenates paired reads into a single fastq file"
    )

    parser.add_argument("-i", "--indir",
                        required=True)
    parser.add_argument("-o", "--outdir",
                        required=True)
    parser.add_argument("-s", "--sbatch",
                        required=True)

    parsed_args = parser.parse_args()
    return parsed_args

def get_files(dir):
    """
    This is expecting files in the format of
    sampleid.R1.fq.gz
    sampleid.R2.fq.gz
    sampleid.unpaired.R1.fq.gz - to ignore
    sampleid.unpaired.R2.fq.gz - to ignore
    <others> will be ignored

    :param dir: directory to search files
    :return: list structured as
            [(sample id, forward file, reverse file)]
            for all sample ids found in the directory
    """
    files = os.listdir(dir)
    fastqs = [file for file in files if "fq" in file]
    pairs = [file for file in fastqs if "unpaired" not in file]
    sample_ids = np.unique([file.split(".")[0] for file in fastqs])

    return_list = []
    for id in sample_ids:
        r1 = pairs[np.where(("R1" in pairs) and (id in pairs))]
        r2 = pairs[np.where(("R1" in pairs) and (id in pairs))]

        return_list.append((id, r1, r2))

    return return_list


def dispatch_zcat(id_r1_r2_list, indir, outdir):
    for id, r1, r2 in id_r1_r2_list:
        in_r1 = os.path.join(indir, r1)
        in_r2 = os.path.join(indir, r2)
        outfile = os.path.join(outdir, f"{id}.concat.fq.gz")

        command = f"zcat {in_r1} {in_r2} | gzip > {outfile}"
        print(f"Executing {command}")
        comp = subprocess.run(command,
                            shell=True)
        print(comp)

if __name__ == "__main__":
    print("Running")
    args = get_args()
    id_r1_r2_list = get_files(args.indir)
    dispatch_zcat(id_r1_r2_list, args.indir, args.outdir)
