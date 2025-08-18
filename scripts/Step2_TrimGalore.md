# Step2: Adapter Trimming and Quality Filtering (Trim Galore!)

This step performs adapter trimming and quality filtering of paired-end FASTQ files using Trim Galore!.
Note: Replace <adapter_sequence> with the actual adapter sequence used in your library preparation kit.
Note: Copy and paste the following commands into your terminal to run them.

## Commands
```
# Run Trim Galore! for adapter trimming and quality filtering
trim_galore --paired \                 # specify paired-end reads
            --quality 20 \             # trim bases with Phred score < 20
            --length 30 \              # discard reads shorter than 30 bases
            --adapter <adapter_sequence> \  # replace with actual adapter sequence
            *_R1.fastq.gz *_R2.fastq.gz     # input FASTQ files (paired-end)
```
## Example (Illumina adapter, 8 threads)
```
trim_galore --paired \
            --quality 20 \
            --length 30 \
            --adapter AGATCGGAAGAGC \
            *_R1.fastq.gz *_R2.fastq.gz
```
