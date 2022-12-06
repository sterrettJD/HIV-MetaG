#!/bin/bash
#SBATCH -p short # Partition or queue.
#SBATCH --job-name=multiqc # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jost9358@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=1 # Run on a single CPU
#SBATCH --mem=4gb # Memory limit
#SBATCH --time=0:15:00 # Time limit hrs:min:sec
#SBATCH --output=/Users/jost9358/MetaG_HIV/slurm_outs/multiqc_%j.out # Standard output and error log
#SBATCH --error=/Users/jost9358/MetaG_HIV/slurm_outs/multiqc_%j.err # %j inserts job number


# README
# PASS TWO PARAMS TO THIS SCRIPT
# -i: your input directory (no slash after)
# -o: your output directory


while getopts i:o: flag
do
    case "${flag}" in
        i) inpath=${OPTARG};;
        o) outpath=${OPTARG};;
    esac
done

source activate /Users/jost9358/miniconda3/envs/multiqc

multiqc $inpath/*_fastqc.zip \
-o $outpath
