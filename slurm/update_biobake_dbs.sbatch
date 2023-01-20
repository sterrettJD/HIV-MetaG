#!/bin/bash
#SBATCH -p short # Partition or queue.
#SBATCH --job-name=db_update # Job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jost9358@colorado.edu
#SBATCH --nodes=1 # Only use a single node
#SBATCH --ntasks=1 # Run on 8 CPUs
#SBATCH --mem=100gb # Memory limit
#SBATCH --time=5:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/update_dbs_%j.out # Standard output and error log
#SBATCH --error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/update_dbs_%j.err # %j inserts job number

source activate humannenv4

humann_databases --download chocophlan full ~/humann_dbs/chocophlan/ --update-config yes
humann_databases --download uniref uniref90_diamond ~/humann_dbs/uniref/ --update-config yes

