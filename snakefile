import pandas as pd

# TODO: document unnecesary f strings
df = pd.read_csv("metadata.txt", sep="\t")
SAMPLES = df["PID"].tolist()

nixing_len = 40
nixing_dir = f"hiv.t32.nix{nixing_len}"



rule all:
# Starting input is data processed through YMP, feels a bit silly to use snakemake to run ymp, which runs snakemake
# Consider the input to be trimmomatic-trimmed forward and reverse reads
    input:
        # Made by nixing short reads
        expand(f"{nixing_dir}/{{sample}}.R1.{nixing_len}.fq.gz", sample=SAMPLES),
        expand(f"{nixing_dir}/{{sample}}.R2.{nixing_len}.fq.gz", sample=SAMPLES),
        # Made by concatenating reads
        # Necessary for steps like humann that require forward and reverse reads to be concatenated
        expand(f"hiv.t32.concat/{{sample}}.concat.fq.gz", sample=SAMPLES),
        # Made by concatenating the nixed reads
        expand(f"hiv.t32.concat.n40/{{sample}}.concat.fq.gz", sample=SAMPLES),
        # Made by nonpareil (WHICH SEEMS TO BE SILENTLY ERRORING OUT MID-RUN)
        expand(f"hiv.t32.concat.n40.nonpareil/{{sample}}.npl", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40.nonpareil/{{sample}}.npo", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40.nonpareil/{{sample}}.npa", sample=SAMPLES),
        # Made by nonpareil when allocated more resources (WHICH SEEMS TO BE SILENTLY ERRORING OUT MID-RUN)
        expand(f"hiv.t32.concat.n40.nonpareil.bigmem/{{sample}}.npl", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40.nonpareil.bigmem/{{sample}}.npo", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40.nonpareil.bigmem/{{sample}}.npa", sample=SAMPLES),

        # Made by Humann (and aggregated by aggregate_humann_outs)
        "hiv.t32.concat.humann/all_pathabundance.tsv",
        "hiv.t32.concat.humann/all_pathcoverage.tsv",
        "hiv.t32.concat.humann/all_genefamilies.tsv",
        "hiv.t32.concat.humann/all_genefamilies_grouped.tsv",
        "hiv.t32.concat.humann/all_genefamilies_grouped.tsv",
        "hiv.t32.concat.humann/all_bugs_list.tsv",
        expand("hiv.t32.concat.humann/{sample}/{sample}.concat_humann_temp/{sample}.concat_metaphlan_bugs_list_v3.tsv",
               sample=SAMPLES),

        # Made by metaspades
        expand(f"hiv.t32.n40.metaspades/{{sample}}/corrected/{{sample}}.R1.{nixing_len}.fq.00.0_0.cor.fastq.gz", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/corrected/{{sample}}.R2.{nixing_len}.fq.00.0_0.cor.fastq.gz", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/corrected/{{sample}}.R_unpaired.{nixing_len}.fq.00.0_0.cor.fastq.gz", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/contigs.fasta", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.fasta", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/contigs.paths", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.paths", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/assembly_graph.fastg", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/assembly_graph_with_scaffolds.gfa", sample=SAMPLES)



rule nix_shortreads:
  input:
      FORWARD=f"hiv.trim_trimmomaticT32/{{sample}}.R1.fq.gz",
      REVERSE=f"hiv.trim_trimmomaticT32/{{sample}}.R2.fq.gz"
  output:
      FORWARD=f"hiv.t32.nix40/{{sample}}.R1.{nixing_len}.fq.gz",
      REVERSE=f"hiv.t32.nix40/{{sample}}.R2.{nixing_len}.fq.gz"
  resources:
        partition="short",
        mem_mb=25000, # MB
        runtime=60 # min
  threads: 8
  run:
      shell("mkdir -p %s" % nixing_dir) #in case this directory doesn't exist. if it does, nothing will be done
      shell(f"bash slurm/nix_shortreads.sh {{input.FORWARD}} {{input.REVERSE}} {nixing_len} {{output.FORWARD}} {{output.REVERSE}}")


rule concat_files:
  input:
      FORWARD=f"hiv.trim_trimmomaticT32/{{sample}}.R1.fq.gz",
      REVERSE=f"hiv.trim_trimmomaticT32/{{sample}}.R2.fq.gz"
  output:
      f"hiv.t32.concat/{{sample}}.concat.fq.gz"
  resources:
        partition="short",
        mem_mb=30000, # MB
        runtime=int(60*2.5) # min
  threads: 1
  run:
      shell("mkdir -p hiv.t32.concat") #in case this directory doesn't exist. if it does, nothing will be done
      shell(f"bash slurm/concat_files.sh -f {{input.FORWARD}} -r {{input.REVERSE}} -o {{output}}")


rule concat_nixed_files:
  input:
      FORWARD=f"hiv.t32.nix40/{{sample}}.R1.{nixing_len}.fq.gz",
      REVERSE=f"hiv.t32.nix40/{{sample}}.R2.{nixing_len}.fq.gz"
  output:
      f"hiv.t32.concat.n40/{{sample}}.concat.fq.gz"
  resources:
        partition="short",
        mem_mb=30000, # MB
        runtime=int(60*2.5) # min
  threads: 1
  run:
      shell("mkdir -p hiv.t32.concat.n40") #in case this directory doesn't exist. if it does, nothing will be done
      shell(f"bash slurm/concat_files.sh -f {{input.FORWARD}} -r {{input.REVERSE}} -o {{output}}")


rule run_nonpareil:
    input:
        "hiv.t32.concat.n40/{sample}.concat.fq.gz"
    output:
        f"hiv.t32.concat.n40.nonpareil/{{sample}}.npl",
        f"hiv.t32.concat.n40.nonpareil/{{sample}}.npo",
        f"hiv.t32.concat.n40.nonpareil/{{sample}}.npa"
    resources:
        partition="short",
        mem_mb=30000, # MB
        runtime=int(60*3) # min
    threads: 16
    conda: "conda_envs/nonpareil.yaml"
    shell:
        """
        mkdir -p hiv.t32.concat.n40.nonpareil
        bash slurm/run_nonpareil.sh -i {input} -o hiv.t32.concat.n40.nonpareil/{wildcards.sample}
        """

# Should probably delete one of these once I figure out the issues with nonpareil
rule run_nonpareil_bigmem:
    input:
        "hiv.t32.concat.n40/{sample}.concat.fq.gz"
    output:
        f"hiv.t32.concat.n40.nonpareil.bigmem/{{sample}}.npl",
        f"hiv.t32.concat.n40.nonpareil.bigmem/{{sample}}.npo",
        f"hiv.t32.concat.n40.nonpareil.bigmem/{{sample}}.npa"
    resources:
        partition="short",
        mem_mb=60000, # MB
        runtime=int(60*5) # min
    threads: 16
    conda: "conda_envs/nonpareil.yaml"
    shell:
        """
        mkdir -p hiv.t32.concat.n40.nonpareil.bigmem
        bash slurm/run_nonpareil.sh -i {input} -o hiv.t32.concat.n40.nonpareil.bigmem/{wildcards.sample}
        """

rule get_biobakery_dbs:
    output:
        "/Users/jost9358/humann_dbs/chocophlan/"
        "/Users/jost9358/humann_dbs/uniref/"
    resources:
        partition="short",
        mem_mb=100000, # MB
        runtime=int(60*5) # min
    # This can also be run using slurm/update_biobake_dbs.sbatch
    threads: 1
    conda: "conda_envs/humann.yaml"
    shell:
        """
        humann_databases --download chocophlan full ~/humann_dbs/chocophlan/ --update-config yes
        humann_databases --download uniref uniref90_diamond ~/humann_dbs/uniref/ --update-config yes
        """

rule run_humann:
  input:
      CHOCO_DB="/Users/jost9358/humann_dbs/chocophlan/",
      UNIREF_DB="/Users/jost9358/humann_dbs/uniref/",
      CONCAT_FILES=f"hiv.t32.concat/{{sample}}.concat.fq.gz"
  output:
      PATHABUND=expand("hiv.t32.concat.humann/{sample}/pathabundance.tsv",
                sample=SAMPLES),
      PATHCOV=expand("hiv.t32.concat.humann/{sample}/pathcoverage.tsv",
                sample=SAMPLES),
      GENEFAMS=expand("hiv.t32.concat.humann/{sample}/genefamilies.tsv",
                sample=SAMPLES),
      BUGSLIST=expand("hiv.t32.concat.humann/{sample}/{sample}.concat_humann_temp/{sample}.concat_metaphlan_bugs_list.tsv",
                sample=SAMPLES)
  resources:
      partition="long",
      mem_mb=int(150*1000), # MB, or 150 GB
      runtime=int(48*60) # min, or 48 hours
  threads: 16
  conda: "conda_envs/humann.yaml"
  shell:
      """
      mkdir -p hiv.t32.concat.humann
      humann -i {input.CONCAT_FILES} -o hiv.t32.concat.humann/{wildcards.sample} --threads 16 --search-mode uniref90
      """


rule aggregate_humann_outs:
    input:
        PATHABUND=expand("hiv.t32.concat.humann/{sample}/pathabundance.tsv",
                sample=SAMPLES),
        PATHCOV=expand("hiv.t32.concat.humann/{sample}/pathcoverage.tsv",
                sample=SAMPLES),
        GENEFAMS=expand("hiv.t32.concat.humann/{sample}/genefamilies.tsv",
                sample=SAMPLES),
        BUGSLIST=expand("hiv.t32.concat.humann/{sample}/{sample}.concat_humann_temp/{sample}.concat_metaphlan_bugs_list.tsv",
                sample=SAMPLES)
    output:
        PATHABUND="hiv.t32.concat.humann/all_pathabundance.tsv",
        PATHCOV="hiv.t32.concat.humann/all_pathcoverage.tsv",
        GENEFAMS="hiv.t32.concat.humann/all_genefamilies.tsv",
        GENEFAMS_GROUPED="hiv.t32.concat.humann/all_genefamilies_grouped.tsv",
        GENEFAMS_GROUPED_NAMED="hiv.t32.concat.humann/all_genefamilies_grouped_named.tsv",
        BUGSLIST="hiv.t32.concat.humann/all_bugs_list.tsv",
        V3_NOAGG_BUGS=expand("hiv.t32.concat.humann/{sample}/{sample}.concat_humann_temp/{sample}.concat_metaphlan_bugs_list_v3.tsv",
                             sample=SAMPLES)

    resources:
        partition="short",
        mem_mb=int(10*1000), # MB, or 10 GB
        runtime=60 # min
    threads: 1
    conda: "conda_envs/humann.yaml"
    shell:
        """
        humann_join_tables -i hiv.t32.concat.humann -o {output.PATHABUND} --file_name pathabundance.tsv --search-subdirectories

        humann_join_tables -i hiv.t32.concat.humann -o {output.PATHCOV} --file_name pathcoverage.tsv --search-subdirectories

        humann_join_tables -i hiv.t32.concat.humann -o {output.GENEFAMS} --file_name genefamilies.tsv --search-subdirectories
        humann_regroup_table -i {output.GENEFAMS} -g uniref90_rxn -o {output.GENEFAMS_GROUPED}
        humann_rename_table -i {output.GENEFAMS_GROUPED} -n metacyc-rxn -o {output.GENEFAMS_GROUPED_NAMED}

        python utils/aggregate_metaphlan_bugslists.py -i hiv.t32.concat.humann -o {output.BUGSLIST}

        python utils/convert_mphlan_v4_to_v3.py -i hiv.t32.concat.humann
        """


rule assemble_metaspades:
    input:
        FORWARD=f"hiv.t32.nix40/{{sample}}.R1.{nixing_len}.fq.gz",
        REVERSE=f"hiv.t32.nix40/{{sample}}.R2.{nixing_len}.fq.gz"
    output:
        CORRECT_R1=f"hiv.t32.n40.metaspades/{{sample}}/corrected/{{sample}}.R1.{nixing_len}.fq.00.0_0.cor.fastq.gz",
        CORRECT_R2=f"hiv.t32.n40.metaspades/{{sample}}/corrected/{{sample}}.R2.{nixing_len}.fq.00.0_0.cor.fastq.gz",
        CORRECT_UNPAIRED=f"hiv.t32.n40.metaspades/{{sample}}/corrected/{{sample}}.R_unpaired.00.0_0.cor.fastq.gz",
        CONTIGS=f"hiv.t32.n40.metaspades/{{sample}}/contigs.fasta",
        SCAFFOLDS=f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.fasta",
        CONTIG_PATHS=f"hiv.t32.n40.metaspades/{{sample}}/contigs.paths",
        SCAFFOLD_PATHS=f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.paths",
        GRAPH=f"hiv.t32.n40.metaspades/{{sample}}/assembly_graph.fastg",
        GRAPH_SCAFFOLDS=f"hiv.t32.n40.metaspades/{{sample}}/assembly_graph_with_scaffolds.gfa"
    resources:
      partition="short",
      mem_mb=int(125*1000), # MB, or 125 GB
      runtime=int(23*60) # min, or 23 hours
    threads: 32
    conda: "conda_envs/metaspades.yaml"
    shell:
        """
        mkdir -p hiv.t32.n40.metaspades
        metaspades.py -o hiv.t32.n40.metaspades/{wildcards.sample} --pe1-1 {input.FORWARD} --pe1-2 {input.REVERSE} --threads 32
        """

# Add in Seqtk for subsampling?

#
# rule map_fastq_to_contigs:
"""
bowtie2-build final.contigs.fa final.contigs
bowtie2 -x final.contigs -1 tara_reads_R1.fastq.gz -2 tara_reads_R2.fastq.gz | \
    samtools view -bS -o tara_to_sort.bam
samtools sort tara_to_sort.bam -o tara.bam
samtools index tara.bam
"""

# MetaBAT2 for binning
#rule run_metabat2:
#    input:
#        SCAFFOLDS=f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.fasta",
#        SORTED_BAM= <NO INPUT YET>
#
#    output:
#        ????? https://nf-co.re/mag/1.1.1/output#spades
#    resources:
#      partition="short",
#      mem_mb=???, # MB, or 125 GB
#      runtime=??? # min, or 23 hours
#    threads: ???
#    conda: "conda_envs/metabat2.yaml"
#    shell:
#        """
#        runMetaBat.sh {input.SCAFFOLDS} {input.SORTED_BAM}
#        """


# PRODIGAL

# CheckM for assessing MAGs

# CheckV for viral MAGs

# EukRep for eukaryyotic MAGs

# dRep for de-replication

# inStrain for strain-level diversity