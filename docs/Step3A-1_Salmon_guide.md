# Step 3A-1: Alignment-free quantification using Salmon

## Tools
- Salmon (e.g., v1.10.0+)

## Inputs
- Reference transcriptome FASTA (`<transcripts.fa>`)
- Paired-end FASTQ files (`<read1.fastq.gz>`, `<read2.fastq.gz>`)

## Outputs
- Salmon index directory (`<index_name>`)
- Per-sample quantification results in `<output_dir>/quant.sf`
  - Contains transcript-level abundance estimates (TPM, NumReads, EffectiveLength, etc.)

## Notes
- Build the Salmon index once per reference transcriptome.
- Run `salmon quant` separately for each sample.
- `quant.sf` files will be used in the next step (gene-level summarization with tximport).
