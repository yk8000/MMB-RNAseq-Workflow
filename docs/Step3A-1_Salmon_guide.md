# Step 3A. Alignment-free quantification using Salmon — Guide

## Purpose
Quantify transcript abundance from paired-end FASTQ without alignment, then summarize to gene level with tximport.

## Inputs
- Paired-end FASTQ files (`*_R1.fastq.gz`, `*_R2.fastq.gz`)
- Reference transcriptome FASTA used to build the Salmon index
- `tx2gene.tsv` (columns: `transcript_id`, `gene_id`) consistent with the same annotation version
- `samples.tsv` (min columns: `sample_id`, `salmon_dir` → per-sample Salmon output directory)

## Outputs
- Per-sample Salmon outputs containing `quant.sf`
- Gene-level tables in `results/tximport/`:
  - `gene_counts.tsv`, `gene_tpm.tsv`, `gene_length.tsv`

## Key options / parameters
- Library type autodetection (`-l A`)
- Selective alignment (`--validateMappings`)
- Threads (`-p`) to match your hardware

## Notes
- Ensure annotation versions match between the index FASTA and `tx2gene.tsv`.
- Run quantification per sample and record output paths in `samples.tsv`.
