#!/bin/bash

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

# make directory for the temp file and output
if [ -d $outpath ]
then
    echo "Directory exists."
else
    echo "Directory $outpath does not exist. Making $outpath now."
    mkdir $outpath
fi


# Unzip fastq
pigz -dc -p 16 $inpath > ${outpath}/temp_unzipped_input.fq

# fastq is recommended for kmer algorithm, so defaulting to those
nonpareil -s ${outpath}/temp_unzipped_input.fq -T kmer -f fastq -b ${outpath} -t 16 #16 threads

# remove the temp file
rm -r ${outpath}
