#!/bin/bash

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
