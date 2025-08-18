# Step 2: Adapter Trimming and Quality Filtering with Trim Galore!

## Tools
- Trim Galore! (e.g., v0.6.7)

## Inputs
- Raw paired-end FASTQ files (`*_R1.fastq.gz`, `*_R2.fastq.gz`)

## Outputs
- Adapter- and quality-trimmed FASTQ files (`*_val_1.fq.gz`, `*_val_2.fq.gz`)

## Notes
- `--quality 20` trims bases with Phred score < 20.
- `--length 30` discards reads shorter than 30 bp after trimming.
- `--adapter <adapter_sequence>` specifies the adapter sequence from your library kit (e.g., Illumina TruSeq: `AGATCGGAAGAGC`).
- Always confirm the adapter sequence in the library preparation manual.
