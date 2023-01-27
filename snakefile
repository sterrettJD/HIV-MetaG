rule all:
# Starting input is data processed through YMP, feels a bit silly to use snakemake to run ymp, which runs snakemake
    input:
        "hiv.trim_trimmomaticT32/{sample}.R1.fq.gz"
        "hiv.trim_trimmomaticT32/{sample}.R2.fq.gz"

# rule nix_shortreads:
#   input:
#       FORWARD = "hiv.trim_trimmomaticT32/{sample}.R1.fq.gz"
#       REVERSE = "hiv.trim_trimmomaticT32/{sample}.R2.fq.gz"
#   output:
#       - nixed forward reads
#       - nixed reverse reads
#   shell:
#       "sbatch nix_shortreads {input.FORWARD} {input.REVERSE} LENGTHPARAM"


# rule concat_files:
#   input:
#       - forward reads
#       - reverse reads
#   output:
#       - concatenated reads
#   shell:
#       "sbatch concat_files -f FORWARD -r REVERSE -o CONCATENATED"

