# Step1: Quality Assessment (FastQC & MultiQC)

This step performs quality assessment of raw FASTQ files.  
**Note**: Copy and paste the following commands into your terminal to run them.  

## Commands

```bash
# Run FastQC for quality assessment
fastqc -t <threads> \        # number of CPU threads to use for parallel processing
       -o <output_dir> \     # directory where FastQC reports will be saved
       --nogroup \           # disable base grouping in per-base plots (useful for long reads)
       -f fastq *.fastq.gz   # input format FASTQ; all gzipped FASTQ files in the current directory

# Aggregate FastQC reports with MultiQC
multiqc <input_dir> -o <output_dir>   # aggregate all FastQC reports from <input_dir> into a summary report
---
# --------------------
# Example (8 threads)
# --------------------
'''bash
       fastqc -t 8
       -o results/qc/fastqc
       --nogroup
       -f fastq data/*.fastq.gz

       multiqc results/qc/fastqc -o results/qc/multiqc
'''
