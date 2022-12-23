# Prevotella rich HIV metagenomics

Sequencing performed end of November 2022

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
   - They took a minimum of ~4 hours and a maximum of 6 hours to run.
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
#### Gene family tables
- **Join:** Used the following command to join the gene family tables `humann_join_tables -i . -o all_genefamilies.tsv --file_name genefamilies.tsv --search-subdirectories`
- **MetaCyc RXN naming:**
  - Used the following command to group the genes into reaction `humann_regroup_table -i all_genefamilies.tsv -g uniref90_rxn -o all_genefamilies_grouped.tsv`
    - This outputs: "Original Feature Count: 19948; Grouped 1+ times: 1016 (5.1%); Grouped 2+ times: 290 (1.5%)"
  - Used the following command to rename the grouped genes `humann_rename_table -i all_genefamilies_grouped.tsv -n metacyc-rxn -o all_genefamilies_grouped_named.tsv`
- **UniRef90 naming:**
  - The UniRef90 mapping file is too small to download with conda/pip, so you need to download it with the following command: `humann_databases --download utility_mapping full $DIR` where `$DIR` is where you want to place the files. I used `../humann_utility/` mapping as my path
    - Read to do that [here](https://forum.biobakery.org/t/rename-table-doesnt-work-with-uniref90/813).
  - 
