#!/bin/bash
#SBATCH -p short # Partition or queue.
#SBATCH --job-name=humann # Job name
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jost9358@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=8 # Run on 8 CPUs
#SBATCH --mem=80gb # Memory limit
#SBATCH --time=23:59:00 # Time limit hrs:min:sec
#SBATCH --output=/Users/jost9358/MetaG_HIV/slurm_outs/humann_%j.out # Standard output and error log
#SBATCH --error=/Users/jost9358/MetaG_HIV/slurm_outs/humann_%j.err # %j inserts job number


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

source activate humannenv

echo $inpath
echo $outpath

humann -i $inpath -o $outpath --threads 8 --search-mode uniref90
