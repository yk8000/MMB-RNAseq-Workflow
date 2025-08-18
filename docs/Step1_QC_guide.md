# Step1: Quality Assessment (FastQC & MultiQC)

## Tools
- FastQC (e.g., v0.11.9)
- MultiQC (e.g., v1.11+)

## Inputs
- Raw FASTQ files (`*.fastq.gz`) located under `--input_dir`

## Outputs
- `results/qc/fastqc/` : per-sample FastQC HTML/ZIP reports  
- `results/qc/multiqc/multiqc_report.html` : aggregated summary
