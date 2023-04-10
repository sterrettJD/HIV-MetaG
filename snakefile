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

        # Made by Humann (and aggregated by aggregate_humann_outs)
        "hiv.t32.concat.humann/all_pathabundance.tsv",
        "hiv.t32.concat.humann/all_pathcoverage.tsv",
        "hiv.t32.concat.humann/all_genefamilies.tsv",
        "hiv.t32.concat.humann/all_genefamilies_grouped.tsv",
        "hiv.t32.concat.humann/all_genefamilies_grouped.tsv",
        "hiv.t32.concat.humann/all_bugs_list.tsv",
        expand("hiv.t32.concat.humann/{sample}/{sample}.concat_humann_temp/{sample}.concat_metaphlan_bugs_list_v3.tsv",
               sample=SAMPLES),

        # made by panphlan (download genomes)
        "prevotella_genomes/Prevotella_copri/Prevotella_copri_pangenome.tsv",
        "prevotella_genomes/Prevotella_copri/panphlan_Prevotella_copri_annot.tsv",
        "prevotella_genomes/Prevotella_copri/Prevotella_copri_pangenome_contigs.fna",
        multiext("prevotella_genomes/Prevotella_copri/Prevotella_copri",
                 ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                 ".rev.1.bt2", ".rev.2.bt2"),
        # made by panphlan (map)
        expand("hiv.t32.p_copri_panphlan/{sample}_p_copri.csv", sample=SAMPLES),

        # Made by metaspades
        expand(f"hiv.t32.n40.metaspades/{{sample}}/corrected/{{sample}}.R1.{nixing_len}.fq.00.0_0.cor.fastq.gz", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/corrected/{{sample}}.R2.{nixing_len}.fq.00.0_0.cor.fastq.gz", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/corrected/{{sample}}.R_unpaired.00.0_0.cor.fastq.gz", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/contigs.fasta", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.fasta", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/contigs.paths", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.paths", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/assembly_graph.fastg", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/assembly_graph_with_scaffolds.gfa", sample=SAMPLES),

        # Made by metaQUAST
#        f"hiv.t32.n40.metaspades.metaQUAST/report.html", # on scaffolds
#        f"hiv.t32.n40.metaspades.metaQUASTc/report.html", # on contigs

        # Made by bowtie2 build
        expand(f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.index{{ext}}", sample=SAMPLES,
               ext=[".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]),
        expand(f"hiv.t32.n40.metaspades/{{sample}}/contigs.index{{ext}}", sample=SAMPLES,
               ext=[".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]),

        # Made by bowtie2 mapping to scaffolds
        expand(f"hiv.t32.n40.metaspades.mapped/{{sample}}_unsorted.bam", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades.mapped/{{sample}}.bam", sample=SAMPLES),

        # Made by bowtie2 mapping to contigs
        expand(f"hiv.t32.n40.metaspades.mappedc/{{sample}}_unsorted.bam", sample=SAMPLES),
        expand(f"hiv.t32.n40.metaspades.mappedc/{{sample}}.bam", sample=SAMPLES),

        # Done file made by metabat2
        expand(f"hiv.t32.n40.metaspades.metabat2/{{sample}}.metabat2done", sample=SAMPLES),

        # CheckM database
        "checkM_db/taxon_marker_sets.tsv",

        # CheckM setup
        expand("hiv.t32.n40.metaspades.metabat2/checkM.copied/{sample}.done",
               sample=SAMPLES),

        # CheckM output
        "hiv.t32.n40.metaspades.metabat2.checkm/checkM.stats.tsv", #scaffolds
        "hiv.t32.n40.metaspades.metabat2.checkmc/checkM.stats.tsv", # contigs

        # CheckV quality summary
        expand(f"hiv.t32.n40.metaspades.checkV/{{sample}}/quality_summary.tsv",
                sample=SAMPLES),

        # Dereplicated genomes from dRep
        "hiv.t32.n40.metaspades.metabat2.checkm.drep/dereplicated_genomes",

        # Classified MAGs from PhyloPhlAn
        "hiv.t32.n40.metaspades.metabat2.checkm.drep.phylophlan_dists",
        "hiv.t32.n40.metaspades.metabat2.checkm.drep.phylophlan_sketches",
        "hiv.t32.n40.metaspades.metabat2.checkm.drep.phylophlan.tsv"


############# PROCESS #############

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


############# ASSESS COVERAGE #############

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

############# RUN BIOBAKERY HUMANN PIPELINE #############

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

############# MAP TO PREVOTELLA COPRI PANGENOME #############
rule download_prevotella_genomes:
    output:
        TSV="prevotella_genomes/Prevotella_copri/Prevotella_copri_pangenome.tsv",
        ANNOT="prevotella_genomes/Prevotella_copri/panphlan_Prevotella_copri_annot.tsv",
        CONTIGS="prevotella_genomes/Prevotella_copri/Prevotella_copri_pangenome_contigs.fna",
        INDEX=multiext(f"prevotella_genomes/Prevotella_copri/Prevotella_copri",
                       ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                       ".rev.1.bt2", ".rev.2.bt2")
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB, or 4 GB
        runtime=int(0.5*60) # min, or 0.5 hours
    threads: 1
    conda:
        "conda_envs/panphlan.yaml"
    shell:
        """
        panphlan_download_pangenome.py -i Prevotella_copri -o prevotella_genomes
        """

rule clean_p_copri_pangenome:
    input:
        ANNOT="prevotella_genomes/Prevotella_copri/panphlan_Prevotella_copri_annot.tsv",
        CONTIGS="prevotella_genomes/Prevotella_copri/Prevotella_copri_pangenome_contigs.fna",
        INDEX=multiext(f"prevotella_genomes/Prevotella_copri/Prevotella_copri",
                       ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                       ".rev.1.bt2", ".rev.2.bt2")
    output:
        ANNOT="prevotella_genomes/Prevotella_copri/panphlan_Prevotella_copri_annot.tsv",
        CONTIGS="prevotella_genomes/Prevotella_copri/Prevotella_copri_pangenome_contigs.fna",
        INDEX=multiext(f"prevotella_genomes/Prevotella_copri/Prevotella_copri",
                       ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                       ".rev.1.bt2", ".rev.2.bt2")
        DONE="prevotella_genomes/Prevotella_copri/genome_cleaned.done"
    resources:
        partition="short",
        mem_mb=int(32*1000), # MB, or 32 GB
        runtime=int(12*60) # min, or 12 hours
    threads: 1
    conda:
        "conda_envs/panphlan.yaml"
    shell:
        """
        panphlan_clean_pangenome.py --species Prevotella_copri \
        --pangenome prevotella_genomes/Prevotella_copri/

        touch {output.DONE}
        """


rule map_panphlan_p_copri:
    input:
        REF_CLEANED="prevotella_genomes/Prevotella_copri/genome_cleaned.done"
        TSV="prevotella_genomes/Prevotella_copri/Prevotella_copri_pangenome.tsv",
        INDEX=multiext(f"prevotella_genomes/Prevotella_copri/Prevotella_copri",
                       ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                       ".rev.1.bt2", ".rev.2.bt2"),
        CONCAT_FILES="hiv.t32.concat/{sample}.concat.fq.gz"
    output:
        MAPPED_CSV="hiv.t32.p_copri_panphlan/{sample}_p_copri.csv"
    resources:
        partition="short",
        mem_mb=int(64*1000), # MB, or 64 GB   TODO: scale down as needed
        runtime=int(23*60) # min, or 23 hours TODO: scale down as needed
    threads: 16
    conda:
        "conda_envs/panphlan.yaml"
    shell:
        """
        mkdir -p hiv.t32.p_copri_panphlan/
        panphlan_map.py -p {input.TSV} \
                    --indexes prevotella_genomes/Prevotella_copri/Prevotella_copri \
                    -i {input.CONCAT_FILES} \
                    --nproc 16 --sam_memory 32 \
                    -o {output.MAPPED_CSV}
        """

# Panphlan profiling

############# ASSEMBLE MAGS #############

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
        mem_mb=int(100*1000), # MB, or 100 GB
        runtime=int(23*60) # min, or 23 hours
    threads: 32
    conda: "conda_envs/metaspades.yaml"
    shell:
        """
        mkdir -p hiv.t32.n40.metaspades
        metaspades.py -o hiv.t32.n40.metaspades/{wildcards.sample} --pe1-1 {input.FORWARD} --pe1-2 {input.REVERSE} --threads 32
        """

############# ASSESS MAGS #############
# QUAST
rule MetaQUAST_scaffolds:
    input:
        SCAFFOLDS=expand(f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.fasta",
                         sample=SAMPLES)
    output:
        REP_HTML=f"hiv.t32.n40.metaspades.metaQUAST/report.html" # HTML version of the report with interactive plots inside
        #https://github.com/ablab/quast/discussions/166
    resources:
        partition="short",
        mem_mb=int(300*1000), # MB, needs a lot apparently
        runtime=int(23.9*60) # min, or just under 24 hours. Icarus takes a long time.
    threads: 8
    conda: "conda_envs/QUAST.yaml"
    shell:
        """
        metaquast.py -o hiv.t32.n40.metaspades.metaQUAST/ hiv.t32.n40.metaspades/*/scaffolds.fasta -t 8 --no-icarus
        """


# QUAST
rule MetaQUAST_contigs:
    input:
        CONTIGS=expand(f"hiv.t32.n40.metaspades/{{sample}}/contigs.fasta",
                         sample=SAMPLES)
    output:
        REP_HTML=f"hiv.t32.n40.metaspades.metaQUASTc/report.html" # HTML version of the report with interactive plots inside
        #https://github.com/ablab/quast/discussions/166
    resources:
      partition="short",
      mem_mb=int(300*1000), # MB, needs a lot
      runtime=int(23.9*60) # min, or just under 24 hours
    threads: 8
    conda: "conda_envs/QUAST.yaml"
    shell:
        """
        metaquast.py -o hiv.t32.n40.metaspades.metaQUASTc/ hiv.t32.n40.metaspades/*/contigs.fasta -t 8 --no-icarus
        """

############# BINNING MAGS #############
rule build_scaffolds_index:
    input:
        SCAFFOLDS=f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.fasta"
    output:
        INDEX=multiext(f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.index",
                       ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                       ".rev.1.bt2", ".rev.2.bt2")
    resources:
        partition="short",
        mem_mb=int(20*1000), # MB, or 20 GB
        runtime=int(1.5*60) # min, or 8 hours
    threads: 8
    conda: "conda_envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build --threads 8 {input.SCAFFOLDS} hiv.t32.n40.metaspades/{wildcards.sample}/scaffolds.index
        """

rule map_fastq_to_scaffolds:
    input:
        INDEX=multiext(f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.index",
                       ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                       ".rev.1.bt2", ".rev.2.bt2"),
        FORWARD=f"hiv.t32.nix40/{{sample}}.R1.{nixing_len}.fq.gz",
        REVERSE=f"hiv.t32.nix40/{{sample}}.R2.{nixing_len}.fq.gz"
    output:
        BAM=f"hiv.t32.n40.metaspades.mapped/{{sample}}_unsorted.bam",
        SORTED_BAM=f"hiv.t32.n40.metaspades.mapped/{{sample}}.bam"
    resources:
        partition="short",
        mem_mb=int(20*1000), # MB, or 20 GB
        runtime=int(8*60) # min, or 8 hours
    threads: 8
    conda: "conda_envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -x hiv.t32.n40.metaspades/{wildcards.sample}/scaffolds.index -1 {input.FORWARD} -2 {input.REVERSE} -p 8 | \
            samtools view -bS -o {output.BAM}
        samtools sort --threads 8 {output.BAM} -o {output.SORTED_BAM}
        samtools index -@ 8 {output.SORTED_BAM}
        """


rule build_contigs_index:
    input:
        CONTIGS=f"hiv.t32.n40.metaspades/{{sample}}/contigs.fasta"
    output:
        INDEX=multiext(f"hiv.t32.n40.metaspades/{{sample}}/contigs.index",
                      ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                      ".rev.1.bt2", ".rev.2.bt2")
    resources:
        partition="short",
        mem_mb=int(20*1000), # MB, or 20 GB
        runtime=int(1.5*60) # min, or 8 hours
    threads: 8
    conda: "conda_envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build --threads 8 {input.CONTIGS} hiv.t32.n40.metaspades/{wildcards.sample}/contigs.index
        """


rule map_fastq_to_contigs:
    input:
        INDEX=multiext(f"hiv.t32.n40.metaspades/{{sample}}/contigs.index",
                      ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                      ".rev.1.bt2", ".rev.2.bt2"),
        FORWARD=f"hiv.t32.nix40/{{sample}}.R1.{nixing_len}.fq.gz",
        REVERSE=f"hiv.t32.nix40/{{sample}}.R2.{nixing_len}.fq.gz"
    output:
        BAM=f"hiv.t32.n40.metaspades.mappedc/{{sample}}_unsorted.bam",
        SORTED_BAM=f"hiv.t32.n40.metaspades.mappedc/{{sample}}.bam"
    resources:
        partition="short",
        mem_mb=int(20*1000), # MB, or 20 GB
        runtime=int(8*60) # min, or 8 hours
    threads: 8
    conda: "conda_envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -x hiv.t32.n40.metaspades/{wildcards.sample}/contigs.index -1 {input.FORWARD} -2 {input.REVERSE} -p 8 | \
            samtools view -bS -o {output.BAM}
        samtools sort -@ 8 {output.BAM} -o {output.SORTED_BAM}
        samtools index -@ 8 {output.SORTED_BAM}
        """


# MetaBAT2 for binning
rule run_metabat2_scaffolds:
    input:
        SCAFFOLDS=f"hiv.t32.n40.metaspades/{{sample}}/scaffolds.fasta",
        SORTED_BAM=f"hiv.t32.n40.metaspades.mapped/{{sample}}.bam"
    output:
        DONE=f"hiv.t32.n40.metaspades.metabat2/{{sample}}.metabat2done"
    resources:
        partition="short",
        mem_mb=int(10*1000), # MB, or 10 GB
        runtime=int(2*60) # min, or 2 hours
    threads: 1
    conda: "conda_envs/metabat2.yaml"
    shell:
        """
        mkdir -p hiv.t32.n40.metaspades.metabat2/{wildcards.sample}
        # Move to this dir so there aren't any issues with each job making a depth.txt file named the same thing
        cd hiv.t32.n40.metaspades.metabat2/{wildcards.sample}
        runMetaBat.sh -m 1500 ../../{input.SCAFFOLDS} ../../{input.SORTED_BAM}
        mv scaffolds.fasta.metabat-bins* bins
        touch ../{wildcards.sample}.metabat2done
        """

rule run_metabat2_contigs:
    input:
        CONTIGS=f"hiv.t32.n40.metaspades/{{sample}}/contigs.fasta",
        SORTED_BAM=f"hiv.t32.n40.metaspades.mappedc/{{sample}}.bam"
    output:
        DONE=f"hiv.t32.n40.metaspades.metabat2c/{{sample}}.metabat2done"
    resources:
        partition="short",
        mem_mb=int(10*1000), # MB, or 10 GB
        runtime=int(2*60) # min, or 2 hours
    threads: 1
    conda: "conda_envs/metabat2.yaml"
    shell:
        """
        mkdir -p hiv.t32.n40.metaspades.metabat2c/{wildcards.sample}
        # Move to this dir so there aren't any issues with each job making a depth.txt file named the same thing
        cd hiv.t32.n40.metaspades.metabat2c/{wildcards.sample}
        runMetaBat.sh -m 1500 ../../{input.CONTIGS} ../../{input.SORTED_BAM}
        mv contigs.fasta.metabat-bins* bins
        touch ../{wildcards.sample}.metabat2done
        """

############# ASSESS MAGS WITH CHECKM #############
# CheckM for assessing MAGs
rule pull_checkM_db:
    output:
        "checkM_db/taxon_marker_sets.tsv"
    resources:
        partition="short",
        mem_mb=int(10*1000), # MB, or 10 GB
        runtime=int(2*60) # min, or 2 hours
    threads: 1
    conda: "conda_envs/checkM.yaml"
    shell:
        """
        mkdir -p checkM_db
        cd checkM_db
        wget https://zenodo.org/record/7401545/files/checkm_data_2015_01_16.tar.gz
        tar -xvf checkm_data_2015_01_16.tar.gz
        cd ..
        checkm data setRoot checkM_db
        """

rule setup_scaffold_bins_for_checkM:
    input:
        BAT2DONE=f"hiv.t32.n40.metaspades.metabat2/{{sample}}.metabat2done"
    output:
        "hiv.t32.n40.metaspades.metabat2/checkM.copied/{sample}.done"
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB, or 8 GB
        runtime=int(1*60) # min, or 1 hour
    threads: 1
    conda: "conda_envs/checkM.yaml"
    shell:
        """
        mkdir -p hiv.t32.n40.metaspades.metabat2/bins_to_derep
        cd hiv.t32.n40.metaspades.metabat2/{wildcards.sample}/bins/

        # copy the bins over and prepend with the sample name so there aren't any with the same name
        for f in bin.*.fa; do cp -v -- "$f" "../../bins_to_derep/{wildcards.sample}.$f"; done

        cd ../../
        mkdir -p checkM.copied/
        touch checkM.copied/{wildcards.sample}.done
        """

rule setup_contig_bins_for_checkM:
    input:
        BAT2DONE=f"hiv.t32.n40.metaspades.metabat2c/{{sample}}.metabat2done"
    output:
        "hiv.t32.n40.metaspades.metabat2c/checkM.copied/{sample}.done"
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB, or 8 GB
        runtime=int(1*60) # min, or 1 hour
    threads: 1
    conda: "conda_envs/checkM.yaml"
    shell:
        """
        mkdir -p hiv.t32.n40.metaspades.metabat2c/bins_to_derep
        cd hiv.t32.n40.metaspades.metabat2c/{wildcards.sample}/bins/

        # copy the bins over and prepend with the sample name so there aren't any with the same name
        for f in bin.*.fa; do cp -v -- "$f" "../../bins_to_derep/{wildcards.sample}.$f"; done

        cd ../../
        mkdir -p checkM.copied/
        touch checkM.copied/{wildcards.sample}.done
        """


# CheckM
rule checkM_scaffolds:
    input:
        DB="checkM_db/taxon_marker_sets.tsv",
        SETUPDONE=expand("hiv.t32.n40.metaspades.metabat2/checkM.copied/{sample}.done",
                        sample=SAMPLES)
    output:
        "hiv.t32.n40.metaspades.metabat2.checkm/checkM.stats.out"
    resources:
        partition="long",
        mem_mb=int(420*1000), # MB, or 350 GB
        runtime=int(30*60) # min, or 30 hours
    threads: 48
    conda: "conda_envs/checkM.yaml"
    shell:
        """
        # Make sure it has our correct DB
        checkm data setRoot checkM_db
        mkdir -p hiv.t32.n40.metaspades.metabat2.checkm
        checkm lineage_wf -t 40 -x fa hiv.t32.n40.metaspades.metabat2/bins_to_derep/ hiv.t32.n40.metaspades.metabat2.checkm > {output}
        """

rule checkM_contigs:
    input:
        DB="checkM_db/taxon_marker_sets.tsv",
        SETUPDONE=expand("hiv.t32.n40.metaspades.metabat2c/checkM.copied/{sample}.done",
                        sample=SAMPLES)
    output:
        "hiv.t32.n40.metaspades.metabat2.checkmc/checkM.stats.out"
    resources:
        partition="long",
        mem_mb=int(420*1000), # MB, or 420 GB
        runtime=int(30*60) # min, or 30 hrs
    threads: 48
    conda: "conda_envs/checkM.yaml"
    shell:
        """
        # Make sure it has our correct DB
        checkm data setRoot checkM_db
        mkdir -p hiv.t32.n40.metaspades.metabat2.checkmc
        checkm lineage_wf -t 40 -x fa hiv.t32.n40.metaspades.metabat2c/bins_to_derep/ hiv.t32.n40.metaspades.metabat2.checkmc > {output}
        """

rule checkM_clean_out:
    input:
        SCAFF="hiv.t32.n40.metaspades.metabat2.checkm/checkM.stats.out",
        CONT="hiv.t32.n40.metaspades.metabat2.checkmc/checkM.stats.out"
    output:
        SCAFFTSV="hiv.t32.n40.metaspades.metabat2.checkm/checkM.stats.tsv",
        SCAFFCSV="hiv.t32.n40.metaspades.metabat2.checkm/checkM.stats.csv",
        CONTTSV="hiv.t32.n40.metaspades.metabat2.checkmc/checkM.stats.tsv",
        CONTCSV="hiv.t32.n40.metaspades.metabat2.checkmc/checkM.stats.csv"
    resources:
        partition="short",
        mem_mb=int(10*1000), # MB, or 10 GB
        runtime=int(1*60) # min, or 1 hours
    threads: 1
    # No conda env - works in just the existing snakemake env
    shell:
        """
        # Move just the stats to a tsv
        sed -n -e '/---------/,$p' {input.SCAFF} | tail -n +2 > {output.SCAFFTSV}
        # convert that tsv to a csv that dRep can use
        # this also converts some column names
        python utils/CheckM_out_to_csv.py -i {output.SCAFFTSV} -o {output.SCAFFCSV}


        # Move just the stats to a tsv (contigs)
        sed -n -e '/---------/,$p' {input.CONT} | tail -n +2 > {output.CONTTSV}
        # convert that tsv to a csv that dRep can use
        # this also converts some column names
        python utils/CheckM_out_to_csv.py -i {output.CONTTSV} -o {output.CONTCSV}
        """

############# DEREPLICATE MAGS #############
# dRep for de-replication
rule dRep_scaffolds:
    input:
        CSV="hiv.t32.n40.metaspades.metabat2.checkm/checkM.stats.csv",
        BINS="hiv.t32.n40.metaspades.metabat2/bins_to_derep"
    output:
        DREPPED=directory("hiv.t32.n40.metaspades.metabat2.checkm.drep/dereplicated_genomes/")
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB, or 8 GB
        runtime=int(2*60) # min, or 2 hours
    threads: 16
    conda: "conda_envs/dRep.yaml"
    shell:
        """
        mkdir -p hiv.t32.n40.metaspades.metabat2.checkm.drep/
        cd hiv.t32.n40.metaspades.metabat2.checkm.drep/

        find ../hiv.t32.n40.metaspades.metabat2/bins_to_derep/ > bins_list.txt

        dRep dereplicate --genomes ../hiv.t32.n40.metaspades.metabat2/bins_to_derep/*.fa \
        -p 16 --debug \
        --genomeInfo ../hiv.t32.n40.metaspades.metabat2.checkm/checkM.stats.csv \
        ./ # specifies current directory as work directory
        """

############# CLASSIFY DEREPLICATED MAGS #############
rule classify_MAGS_phylophlan:
    input:
        SGB="hiv.t32.n40.metaspades.metabat2.checkm.drep/dereplicated_genomes/"
    output:
        directory("hiv.t32.n40.metaspades.metabat2.checkm.drep.phylophlan_dists"),
        directory("hiv.t32.n40.metaspades.metabat2.checkm.drep.phylophlan_sketches"),
        "hiv.t32.n40.metaspades.metabat2.checkm.drep.phylophlan.tsv"
    resources:
        partition="short",
        mem_mb=int(23*1000), # MB, or 24 GB
        runtime=int(12*60) # min, or 23 hours TODO: see if needs to be scaled up/down
    threads: 1
    conda: "conda_envs/phylophlan.yaml"
    retries: 3
    shell:
        """
        # will pull the July 2020 database from Pasoli et al. 2019 https://doi.org/10.1016/j.cell.2019.01.001
        phylophlan_metagenomic \
            -i {input.SGB} \
            -o {output} \
            --nproc 1 \
            -n 1 \
            -d SGB.Jul20 \
            --verbose
        """

############# ASSESS VIRAL MAGS #############
# CheckV for viral MAGs
rule pull_checkV_db:
    output:
        "checkv-db-v1.5/genome_db/checkv_reps.faa" # just one of the files
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB, or 4 GB
        runtime=int(0.5*60) # min
    threads: 1
    conda: "conda_envs/checkV.yaml"
    shell:
        """
        checkv download_database ./
        export CHECKVDB=checkv-db-v1.5
        """

rule checkV:
    input:
        DB="checkv-db-v1.5/genome_db/checkv_reps.faa",
        CONTIGS=f"hiv.t32.n40.metaspades/{{sample}}/contigs.fasta"
    output:
        MAIN=directory(f"hiv.t32.n40.metaspades.checkV/{{sample}}/"),
        QUAL=f"hiv.t32.n40.metaspades.checkV/{{sample}}/quality_summary.tsv"
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB, or 16 GB
        runtime=int(23*60) # min, or 23 hours
    threads: 16
    conda: "conda_envs/checkV.yaml"
    shell:
        """
        export CHECKVDB=checkv-db-v1.5
        mkdir -p hiv.t32.n40.metaspades.checkV/
        checkv end_to_end {input.CONTIGS} hiv.t32.n40.metaspades.checkV/{wildcards.sample}/ -t 16
        """


# EukRep for eukaryyotic MAGs

# inStrain for strain-level diversity?

# Use https://instrain.readthedocs.io/en/latest/user_manual.html#gene-annotation to look at annotation
