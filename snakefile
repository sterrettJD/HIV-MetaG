import pandas as pd

df = pd.read_csv("metadata.txt", sep="\t")
SAMPLES = df["PID"].tolist()

nixing_len = 40
nixing_dir = f"hiv.t32.nix{nixing_len}"



rule all:
# Starting input is data processed through YMP, feels a bit silly to use snakemake to run ymp, which runs snakemake
# Consider the input to be trimmomatic-trimmed forward and reverse reads
    input:
        expand(f"{nixing_dir}/{{sample}}.R1.{nixing_len}.fq.gz", sample=SAMPLES),
        expand(f"{nixing_dir}/{{sample}}.R2.{nixing_len}.fq.gz", sample=SAMPLES),
        expand(f"hiv.t32.concat/{{sample}}.concat.fq.gz", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40/{{sample}}.concat.fq.gz", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40.nonpariel/{{sample}}.npl", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40.nonpariel/{{sample}}.npo", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40.nonpariel/{{sample}}.npa", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40.nonpariel.bigmem/{{sample}}.npl", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40.nonpariel.bigmem/{{sample}}.npo", sample=SAMPLES),
        expand(f"hiv.t32.concat.n40.nonpariel.bigmem/{{sample}}.npa" sample=SAMPLES)


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
        runtime=60, # min
        tasks=8,
        slurm_extra="--error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.err --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.out --mail-type=END --mail-user=jost9358@colorado.edu"
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
        runtime=int(60*2.5), # min
        tasks=1,
        slurm_extra="--error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.err --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.out --mail-type=END --mail-user=jost9358@colorado.edu"
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
        runtime=int(60*2.5), # min
        tasks=1,
        slurm_extra="--error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.err --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.out --mail-type=END --mail-user=jost9358@colorado.edu"
  run:
      shell("mkdir -p hiv.t32.concat.n40") #in case this directory doesn't exist. if it does, nothing will be done
      shell(f"bash slurm/concat_files.sh -f {{input.FORWARD}} -r {{input.REVERSE}} -o {{output}}")


rule run_nonpariel:
    input:
        "hiv.t32.concat.n40/{sample}.concat.fq.gz"
    output:
        f"hiv.t32.concat.n40.nonpariel/{{sample}}.npl",
        f"hiv.t32.concat.n40.nonpariel/{{sample}}.npo",
        f"hiv.t32.concat.n40.nonpariel/{{sample}}.npa"
    resources:
        partition="short",
        mem_mb=30000, # MB
        runtime=int(60*3), # min
        tasks=16,
        slurm_extra="--error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.err --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.out --mail-type=END --mail-user=jost9358@colorado.edu"
    run:
        shell("mkdir -p hiv.t32.concat.n40.nonpariel")
        shell("bash slurm/run_nonpariel.sh -i {input} -o hiv.t32.concat.n40.nonpariel/{wildcards.sample}")

rule run_nonpariel_bigmem:
    input:
        "hiv.t32.concat.n40/{sample}.concat.fq.gz"
    output:
        f"hiv.t32.concat.n40.nonpariel.bigmem/{{sample}}.npl",
        f"hiv.t32.concat.n40.nonpariel.bigmem/{{sample}}.npo",
        f"hiv.t32.concat.n40.nonpariel.bigmem/{{sample}}.npa"
    resources:
        partition="short",
        mem_mb=60000, # MB
        runtime=int(60*5), # min
        tasks=16,
        slurm_extra="--error=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.err --output=/scratch/Users/jost9358/HIV-MetaG/slurm_outs/nixshort_%j.out --mail-type=END --mail-user=jost9358@colorado.edu"
    run:
        shell("mkdir -p hiv.t32.concat.n40.nonpariel.bigmem")
        shell("bash slurm/run_nonpariel.sh -i {input} -o hiv.t32.concat.n40.nonpariel.bigmem/{wildcards.sample}")


"""
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
"""
