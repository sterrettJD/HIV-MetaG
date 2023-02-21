#!/bin/bash
#SBATCH -p short # Partition or queue.
#SBATCH --job-name=metaQUAST # Job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jost9358@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=8 # Run on a single CPU
#SBATCH --mem=50gb # Memory limit
#SBATCH --time=08:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/metaQUAST_%j.out # Standard output and error log
#SBATCH --error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/metaQUAST_%j.err # %j inserts job number


source activate /Users/jost9358/miniconda3/envs/QUAST

metaquast.py -o quast_test hiv.t32.n40.metaspades/DM001/scaffolds.fasta -t 8