# Prevotella rich HIV/MSM shotgun metagenomics

Sequencing performed end of November 2022

Main sections:
1. Processing/QC with YMP
2. Processing in Snakemake
    - Removing short reads (for NonPareil) TODO: add to readme
    - Concatenating paired reads
    - Running HUMANN pipeline
    - Aggregating HUMANN results
    - Pull in other datasets to assess strain-level differences in Prevotella
    - Assessing metagenome coverage (NonPareil) TODO: add to readme
    - Assembling genomes (MetaSPAdes) TODO: add to readme, perform
    - Assessing strain-level variation (inStrain, other tools too?) TODO: all
    - Assess SNPs in Prevotella genomes (MIDAS) TODO: all
3. Analysis of taxonomic/functional data
4. Analysis of strain-level differences (TODO: all)

## Data processing with YMP
### Creation of YMP metadata
 - `utils/demux_to_ymp_metadata.py` was used to convert the genomics core formatted demux csv to YMP metadata with filepaths
### YMP metadata
 - `ymp.yml` contains the YMP file for processing data
### Initial quality checking
 - All raw data files were run through fastqc using `ymp submit hiv.myqc`, which runs the YMP stage `qc_fastq`
 - QC reports were concatenated via MultiQC using the script `sbatch slurm/run_multiqc.sbatch -i hiv.myqc -o multiqc_report/`
 - Results
   - Most checks passed, except one sample (LG001.R2) failed the per base sequence content, and 6 samples had overrepresented sequences. I'm not too concerned about these, since hopefully the trimming and mapping with get rid of problematic reads (e.g., maybe human reads are overrepresented). 
   - ALL samples failed the adapter content check. Trimmomatic should be able to fix this, as it will remove the adapters.
### Read trimming and subsequent QC
  - Raw reads were trimmed using Trimmomatic through YMP, using `ymp submit hiv.trim_trimmomatic` 
  - Trimmed reads were then quality checked again to ensure the removal of adapters using `ymp submit hiv.trim_trimmomaticT32.qc_fastqc`
  - Trimmed read quality files were concatenated using MultiQC via `sbatch slurm/run_multiqc.sbatch -i hiv.trim_trimmomaticT32.qc_fastqc -o trimmed_multiqc_out/`
  - With trimming using the T32 parameter, all samples passed all checks
    - Many of the reverse fastqs have some moderate quality reads, so maybe that could be something to address at some point

## Data processing outside of YMP
### Snakemake
- All important following steps should be in `snakefile`
- Run using `snakemake --profile slurm`
  - This requires you to have a profile setup called `slurm` in your .config.
  Mine is located at `~/.config/snakemake/slurm/` on my cluster. 
  I used a [Cookiecutter snakemake slurm](https://github.com/Snakemake-Profiles/slurm). 
- If jobs are running when FIJI times out, run `snakemake --unlock` next time you need to run snakemake. This is because snakemake will lock the directory if it closes unexpectedly.
- A dag of jobs can be made using `snakemake --dag | dot -Tsvg > dag.svg`

### Concatenating paired reads for HUMAnN
 - For HUMAnN, both forward and paired reads need to be in a single fastq file.
 - There is no need to merge them, since output abundances are normalized by kilobase.
 - I concatenated using a script `utils/concat_paired_reads.py`, which calls a simple `.sbatch` in the `slurm/` directory to schedule this job with slurm
 - Specifically, I ran `python utils/concat_paired_reads.py -i hiv.trim_trimmomaticT32 -o hiv.t32.concat/ -s slurm/concat_files.sbatch`

### Running HUMAnN
 - HUMAnN is run via calling `utils/run_humann.py`
 - `utils/run_humann.py` dispatches the sbatch script `slurm/run_humann.sbatch`
 - I did it using the following command `python utils/run_humann.py -i hiv.t32.concat/ -o hiv.t32.concat.humann -s slurm/run_humann.sbatch`
 - This executed one slurm job for every sample. Each job was allocated 8 cores, and here are the stats for reference: 
   - They took a minimum of ~20 hours and a maximum of 36 hours to run.
   - The input files were 5-8 GB.
   - The output files totaled to 2 TB.

#### Setting up HUMAnN (notes)
 - Had some issues getting HUMAnN to work, see [here](https://forum.biobakery.org/t/humann-conda-installation-dependency-issues/4557?u=sterrettjd)
 - In short, some advice:
   - Create a new HUMAnN environment, and set conda channel priorities to prioritize biobakery and conda-forge
   - Use mamba, not conda to install HUMAnN. Conda would take over an hour to try to solve the environment
   - Make sure you're using DIAMOND >= 2.0.15. Previous releases of DIAMOND had bugs, and HUMAnN will quit after about 30 min if you don't. I'm not sure why conda would install an old release of DIAMOND, but it also had issues with the DIAMOND dependencies, but mamba was able to fix those.
   - Update MetaPhlan to the latest release of version 4. HUMAnN install was using version 3.something.
   - Run `metaphlan --install` to install databases and build the index for mapping. This takes multiple hours, so plan accordingly.
     - You can also find prebuilt indexes online if you don't want to build them on your own computer.
   
### Concatenating HUMAnN results
- HUMAnN output a directory for each sample, with the gene families, pathway abundances, and pathway coverage for that sample
- There are also temporary files in a temp subdirectory for each sample, which contain the MetaPhlan mapping outputs
- I navigated to `hiv.t32.concat.humann` to run the following commands.
- The concatenated files have been uploaded to GitHub in the `hiv.t32.concat.humann` directory.

#### Gene family tables
- **Join:** Used the following command to join the gene family tables `humann_join_tables -i . -o all_genefamilies.tsv --file_name genefamilies.tsv --search-subdirectories`
- **MetaCyc RXN naming:**
  - Used the following command to group the genes into reaction `humann_regroup_table -i all_genefamilies.tsv -g uniref90_rxn -o all_genefamilies_grouped.tsv`
    - This outputs: "Original Feature Count: 1701936; Grouped 1+ times: 183290 (10.8%); Grouped 2+ times: 51421 (3.0%)"
  - Used the following command to rename the grouped genes `humann_rename_table -i all_genefamilies_grouped.tsv -n metacyc-rxn -o all_genefamilies_grouped_named.tsv`
- **UniRef90 naming:**
  - The UniRef90 mapping files are too big to download with conda/pip, so you need to download it with the following command: `humann_databases --download utility_mapping full $DIR` where `$DIR` is where you want to place the files. I used `../humann_utility/` mapping as my path
    - Read to do that [here](https://forum.biobakery.org/t/rename-table-doesnt-work-with-uniref90/813).
#### Pathway abundance tables
- **Join:** Used the following command to join the gene family tables `humann_join_tables -i . -o all_pathabundance.tsv --file_name pathabundance.tsv --search-subdirectories`

#### Pathway coverage tables
- **Join:** Used the following command to join the gene family tables `humann_join_tables -i . -o all_pathcoverage.tsv --file_name pathcoverage.tsv --search-subdirectories`

#### MetaPhlan taxonomic profiles
- **Join:** Used the following command to join the taxonomic profiles: `python ../utils/aggregate_metaphlan_bugslists.py -i . -o all_bugs_list.tsv`
- This calls a script in utils that I wrote to aggregate the metaphlan relative abundance data into a single table
- I set it up to work with the naming scheme generated by all steps so far, but it could probably modified if desired.

## Analysis of HUMAnN-generated data
### Pavian
#### Formatting
- Pavian currently doesn't work with MetaPhlan v4 formatted outputs, only v3. I've opened an issue [here](https://github.com/fbreitwieser/pavian/issues/99) with a suggested fix, so hopefully this will be updated. 
- The script `utils/convert_mphlan_v4_to_v3.py` reformats metaphlan v4 reports to trick Pavian into thinking they're v3 reports, which gets around the issue for now.
- I navigate to the humann output directory and run it using `python ../utils/convert_mphlan_v4_to_v3.py -i .`
#### Running Pavian
- To run pavian, see `analysis/run_pavian.R`, it's basically just `pavian::runApp(port=5000)`.
- Pavian will open a new window for the shiny app 
- Select "Read from server" and type in the path`/path_to_repo/HIV-MetaG/hiv.t32.concat.humann/*/*/*_v3.tsv`, then click the button to read those in.

### Taxonomy
- `analysis/Taxonomy.Rmd` contains taxa barplots and comparison of 16S Prevotella to metagenome Prevotella relative abundances.

#### Procrustes analysis
- How well do samples match?
- see if any samples match up better from the permutations

#### Greengenes vs SILVA vs new Greengenes on 16S

### Functional profiles

## Estimating metagenome coverage
### Nonpareil
- See documentation [here](https://nonpareil.readthedocs.io/en/latest/)
- Installed via conda on Jan 26, 2023
#### Running Nonpareil
- Run in slurm using `slurm/run_nonpariel.sbatch`
- 

## Assembling metagenome-assembled genomes (MAGs)
- Don't coassemble, our reads are deep anyway, and strain level variation might break that up
- presence/absence

### MetaSPAdes

### Prodigal

## Pangenome
- see [r gnavus paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0490-5#Sec2)

## Strain level (SNPs)
- see [Katherine Pollard SNP paper] (https://genome.cshlp.org/content/26/11/1612.short)
- [MIDAS](https://github.com/snayfach/MIDAS)

### Other datasets

