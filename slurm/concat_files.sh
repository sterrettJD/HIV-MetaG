#!/bin/bash
# README
# PASS THREE PARAMS TO THIS SCRIPT
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
