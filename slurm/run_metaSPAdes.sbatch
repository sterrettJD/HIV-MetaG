#!/bin/bash
#SBATCH -p long # Partition or queue.
#SBATCH --job-name=metaspades # Job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jost9358@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=32 # Run on a single CPU
#SBATCH --mem=250gb # Memory limit
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/metaspades_%j.out # Standard output and error log
#SBATCH --error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/metaspades_%j.err # %j inserts job number

# README
# PASS THREE PARAMS TO THIS SCRIPT
# -f: forward reads
# -r: reverse reads
# -o: your output directory


while getopts f:r:o: flag
do
    case "${flag}" in
        f) forward=${OPTARG};;
        r) reverse=${OPTARG};;
        o) outpath=${OPTARG};;
    esac
done

metaspades.py -o $outpath --pe1-1 $forward --pe1-2 $reverse --threads 32
