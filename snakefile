rule all:
# Starting input is data processed through YMP, feels a bit silly to use snakemake to run ymp, which runs snakemake
# Consider the input to be trimmomatic-trimmed forward and reverse reads
    input:
        "hiv.trim_trimmomaticT32/{sample}.R1.fq.gz"
        "hiv.trim_trimmomaticT32/{sample}.R2.fq.gz"

# TODO: FORMAT CORRECTLY
# rule nix_shortreads:
#   input:
#       FORWARD = "hiv.trim_trimmomaticT32/{sample}.R1.fq.gz"
#       REVERSE = "hiv.trim_trimmomaticT32/{sample}.R2.fq.gz"
#   output:
#       - nixed forward reads (formatted as FORWARD.LENGTHPARAM)
#       - nixed reverse reads (REVERSE.LENGTHPARAM)
#   shell:
#       "sbatch slurm/nix_shortreads.sbatch {input.FORWARD} {input.REVERSE} LENGTHPARAM"

# TODO: FORMAT CORRECTLY
# rule concat_files:
#   input:
#       - forward reads
#       - reverse reads
#   output:
#       - concatenated reads
#   shell:
#       "sbatch slurm/concat_files.sbatch -f FORWARD -r REVERSE -o CONCATENATED"

# TODO: FORMAT CORRECTLY
# rule concat_nixed_files:
#   input:
#       - forward reads
#       - reverse reads
#   output:
#       - concatenated reads
#   shell:
#       "sbatch slurm/concat_files.sbatch -f FORWARD -r REVERSE -o CONCATENATED"

rule get_biobakery_dbs:
    output:
        "~/humann_dbs/chocophlan/"
        "~/humann_dbs/uniref/"
    shell:
        "sbatch slurm/update_biobake_dbs.sbatch"

# TODO: FORMAT CORRECTLY
# rule run_humann:
#   input:
#       - concatenated trimmed reads
#   output:
#       - outdir/path abundance
#       - outdir/path coverage
#       - outdir/gene families
#       - outdir/sample_temp/bugslist.tsv
#   shell:
#       "sbatch slurm/run_humann.sbatch -i CONCATENATED -o hiv.t32.concat.humann.full.{sample}"

# TODO: manage conda environment
rule aggregate_bugslists:
    input:
        - all of the bugslists
    output:
        hiv.t32.concat.humann.full/all_bugslist.tsv
    shell:
        "python utils/aggregate_metaphlan_bugslists.py -i hiv.t32.concat.humann.full -o {output}"

# TODO: manage conda environment
rule aggregate_humann_pathcoverage:
    input:
        - all of the humann pathcoverage
    output:
        hiv.t32.concat.humann.full/all_pathcoverage.tsv
    shell:
        "source activate humannenv4"
        "humann_join_tables -i hiv.t32.concat.humann.full -o {output} --file_name pathcoverage.tsv --search-subdirectories"

# TODO: manage conda environment
rule aggregate_humann_pathabundance:
    input:
        - all of the humann pathabundance
    output:
        hiv.t32.concat.humann.full/all_pathabundance.tsv
    shell:
        "source activate humannenv4"
        "humann_join_tables -i hiv.t32.concat.humann.full -o {output} --file_name pathabundance.tsv --search-subdirectories"

# TODO: manage conda environment
rule aggregate_humann_genefams:
    input:
        - all of the humann genefamilies
    output:
        hiv.t32.concat.humann.full/all_genefamilies.tsv
    shell:
        "source activate humannenv4"
        "humann_join_tables -i hiv.t32.concat.humann.full -o {output} --file_name genefamilies.tsv --search-subdirectories"

# TODO: manage conda environment
rule group_humann_genefams:
    input:
        hiv.t32.concat.humann.full/all_genefamilies.tsv
    output:
        hiv.t32.concat.humann.full/all_genefamilies_grouped.tsv
    shell:
        "source activate humannenv4"
        "humann_regroup_table -i {input} -g uniref90_rxn -o {output}"

# TODO: manage conda environment
rule rename_humann_genefams_metacyc:
    input:
        - hiv.t32.concat.humann.full/all_genefamilies_grouped.tsv
    output:
        - hiv.t32.concat.humann.full/all_genefamilies_grouped_named.tsv
    shell:
        "source activate humannenv4"
        "humann_rename_table -i {input} -n metacyc-rxn -o {output}"

# TODO: ADD OTHER NONPARIEL FILES
rule run_nonpariel:
    input:
        - nixed trimmed concatenated reads
    output:
        - hiv.t32.nix.concat.nonpariel/{sample}/{sample}.npl
        - OTHER NONPARIEL FILES
    shell:
        "sbatch slurm/run_nonpariel.sbatch -i {input} -o hiv.t32.nix.concat.nonpariel/{sample}/"

rule assemble_metaspades:
    input:
        - nixed FORWARD
        - nixed REVERSE
    output:
        # TODO: figure out structure of metaSPAdes output
    shell:
        "sbatch slurm/run_metaSPAdes.sbatch -f FORWARD -r REVERSE -o OUTPUT"