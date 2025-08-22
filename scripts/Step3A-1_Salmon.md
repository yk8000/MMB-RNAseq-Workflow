# Step 3A. Alignment-free quantification using Salmon

This step builds a Salmon index and quantifies paired-end reads.  


## Commands

```bash
# Build Salmon index
salmon index \
  -t <transcripts.fa>   # reference transcriptome FASTA
  -i <index_name>       # Salmon index directory name/path
  -p <threads>          # number of threads for indexing

# Quantification (per sample)
salmon quant \
  -i <index_name>       # Salmon index directory
  -l A                  # auto-detect library type
  -1 <read1.fastq.gz>   # paired-end read1 FASTQ
  -2 <read2.fastq.gz>   # paired-end read2 FASTQ
  -p <threads>          # number of threads for quant
  --validateMappings    # enable selective alignment
  -o <output_dir>       # output dir (contains quant.sf)
```

## Example

```
# Build index
salmon index \
  -t ref/transcripts.fa \
  -i ref/salmon_index \
  -p 8

# Quantify one sample
salmon quant \
  -i ref/salmon_index \
  -l A \
  -1 data/S1_R1.fastq.gz \
  -2 data/S1_R2.fastq.gz \
  -p 8 \
  --validateMappings \
  -o results/salmon/S1
```

