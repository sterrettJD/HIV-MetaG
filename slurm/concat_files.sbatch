#!/bin/bash
#SBATCH -p short # Partition or queue.
#SBATCH --job-name=concat # Job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jost9358@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=1 # Run on a single CPU
#SBATCH --mem=30gb # Memory limit
#SBATCH --time=2:30:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/concat_%j.out # Standard output and error log
#SBATCH --error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/concat_%j.err # %j inserts job number


# README
# PASS TWO PARAMS TO THIS SCRIPT
# -f: forward reads
# -r: reverse reads
# -o: your output file


while getopts f:r:o: flag
do
    case "${flag}" in
        f) forward=${OPTARG};;
        r) reverse=${OPTARG};;
        o) outpath=${OPTARG};;
    esac
done

zcat $forward $reverse | gzip > $outpath
