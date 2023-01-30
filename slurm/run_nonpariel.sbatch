#!/bin/bash
#SBATCH -p short # Partition or queue.
#SBATCH --job-name=nonpariel # Job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jost9358@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=16 # Run on a single CPU
#SBATCH --mem=20gb # Memory limit
#SBATCH --time=23:59:59 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nonpariel_%j.out # Standard output and error log
#SBATCH --error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nonpariel_%j.err # %j inserts job number

# README
# PASS TWO PARAMS TO THIS SCRIPT
# -i: your fastq file
# -o: your output directory

while getopts i:o: flag
do
    case "${flag}" in
        i) inpath=${OPTARG};;
        o) outpath=${OPTARG};;
    esac
done

# Activate conda enviroment
source activate nonparielenv

# make directory for the temp file and output
if [ -d $outpath ] 
then
    echo "Directory exists." 
else
    echo "Directory does not exist. Making Outpath now."
    mkdir $outpath
fi


# Unzip fastq
zcat $inpath > ${outpath}/temp_unzipped_input.fq

# fastq is recommended for kmer algorithm, so defaulting to those
nonpareil -s ${outpath}/temp_unzipped_input.fq -T kmer -f fastq -b ${outpath} -t 16 #16 threads

# remove the temp file
rm -r ${outpath}
