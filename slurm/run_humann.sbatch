#!/bin/bash
#SBATCH -p long # Partition or queue.
#SBATCH --job-name=humann # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jost9358@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=16 # Run on 16 CPUs
#SBATCH --mem=150gb # Memory limit
#SBATCH --time=47:59:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/humann_%j.out # Standard output and error log
#SBATCH --error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/humann_%j.err # %j inserts job number


# README
# PASS TWO PARAMS TO THIS SCRIPT
# -i: concatenated reads
# -o: your output directory


while getopts i:o: flag
do
    case "${flag}" in
        i) inpath=${OPTARG};;
        o) outpath=${OPTARG};;
    esac
done

source activate humannenv4

echo $inpath
echo $outpath

humann -i $inpath -o $outpath --threads 16 --search-mode uniref90
