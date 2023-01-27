#!/bin/bash
#SBATCH -p short # Partition or queue.
#SBATCH --job-name=nixshort # Job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jost9358@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=8 # Run on a single CPU
#SBATCH --mem=25gb # Memory limit
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.out # Standard output and error log
#SBATCH --error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.err # %j inserts job number


# This removes reads of a below a certain length from paired read files in fastq format (e.g., R1 and R2 from the same library)
# Graciously adapted from https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/27218-removing-short-reads-from-paired-end-fastqs?t=31845
# Usage: $ bash nixshorts_PE [input fastqR1] [input fastqR2] [minimum read length to keep]

# PROCESS:

#1a. Start with inputs
R1fq=$1
R2fq=$2
minlen=$3

#1b.
pigz $R1fq --decompress -p 8 > temp.unzipped.R1
pigz $R2fq --decompress -p 8 > temp.unzipped.R2

#2. Find all entries with read length less than minimum length and print line numbers, for both R1 and R2
awk -v min=$minlen '{if(NR%4==2) if(length($0)<min) print NR"\n"NR-1"\n"NR+1"\n"NR+2}' temp.unzipped.R1 > temp.lines1
awk -v min=$minlen '{if(NR%4==2) if(length($0)<min) print NR"\n"NR-1"\n"NR+1"\n"NR+2}' temp.unzipped.R2 >> temp.lines1

#3. Combine both line files into one, sort them numerically, and collapse redundant entries
sort -n temp.lines1 | uniq > temp.lines
rm temp.lines1

#4. Remove the line numbers recorded in "lines" from both fastqs
awk 'NR==FNR{l[$0];next;} !(FNR in l)' temp.lines temp.unzipped.R1 > $R1fq.$minlen
awk 'NR==FNR{l[$0];next;} !(FNR in l)' temp.lines temp.unzipped.R2 > $R2fq.$minlen

#5. Remove temp files and zip outputs
rm temp.lines
rm temp.unzipped.R1
rm temp.unzipped.R2
pigz $R1fq.$minlen --processes 8
pigz $R2fq.$minlen --processes 8

#6. Conclude
echo "Pairs shorter than $minlen bases removed from $R1fq and $R2fq"