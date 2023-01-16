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
        prog="Run Humann",
        description="Dispatches sbatch to run Humann for each sample"
    )

    parser.add_argument("-i", "--indir",
                        required=True)
    parser.add_argument("-o", "--outdir",
                        required=True)
    parser.add_argument("-s", "--sbatch",
                        required=True)

    parsed_args = parser.parse_args()
    return parsed_args


def make_outdir(outdir):
    if os.path.isdir(outdir) is False:
        os.mkdir(outdir)

    return None


def get_files(dir):
    """
    This is expecting files in the format of
    sampleid.concat.fq.gz
    <others> will be ignored

    :param dir: directory to search files
    :return: list structured as
            [(sample id, file), (sample id, file)]
            for all sample ids found in the directory
    """
    files = os.listdir(dir)
    fastqs = [file for file in files if "fq" in file]
    sample_ids = np.unique([file.split(".")[0] for file in fastqs])

    return_list = []
    for id in sample_ids:
        reads = [file for file in fastqs if id in file][0]
        return_list.append((id, reads))

    return return_list


def dispatch_humann(id_read_list, indir, outdir, slurm_script):
    for id, read in id_read_list:
        in_read = os.path.join(indir, read)
        outfile = os.path.join(outdir, f"{id}")

        command = f"sbatch {slurm_script} -i {in_read} -o {outfile}"
        print(f"Executing {command}")
        comp = subprocess.run(command,
                              shell=True)
        print(comp)


if __name__ == "__main__":
    args = get_args()
    id_r1_r2_list = get_files(args.indir)
    make_outdir(args.outdir)
    dispatch_humann(id_r1_r2_list, args.indir, args.outdir, args.sbatch)
