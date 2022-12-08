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
