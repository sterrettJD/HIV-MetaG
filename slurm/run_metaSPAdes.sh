#!/bin/bash

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
