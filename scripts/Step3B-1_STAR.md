# Step 3- Option B. Alignment-based quantification using STAR + featureCounts
## 1) Genome Indexing and Read Alignment with STAR

**Note:** Copy & paste the following commands into your terminal and replace placeholders `<>`.

## Commands

```bash
#------------------------------------------------------------------------
# Build STAR genome index
STAR \
  --runThreadN <threads> \                 # number of CPU threads
  --runMode genomeGenerate \               # run STAR in genome indexing mode
  --genomeDir <genome_index_dir> \         # directory to store the genome index
  --genomeFastaFiles <reference_genome.fasta> \  # reference genome FASTA
  --sjdbGTFfile <gene_annotation.gtf>      # gene annotation (GTF) for splice junctions

# Align reads with STAR
STAR \
  --runThreadN <threads> \                 # number of CPU threads
  --genomeDir <genome_index_dir> \         # path to the pre-built STAR index
  --readFilesIn <read1.fastq.gz> <read2.fastq.gz> \  # paired-end FASTQ files
  --outSAMtype BAM SortedByCoordinate \    # output coordinate-sorted BAM
  --outFileNamePrefix <output_prefix>      # prefix for output filenames
#------------------------------------------------------------------------
```

## Example
```
# Build STAR index (GRCh38 example)
STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir ref/STAR_index \
  --genomeFastaFiles ref/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile ref/gencode.v46.annotation.gtf

# Align one sample (gzipped FASTQ)
# If inputs are .fastq.gz, add: --readFilesCommand zcat
STAR \
  --runThreadN 8 \
  --genomeDir ref/STAR_index \
  --readFilesIn data/S1_R1.fastq.gz data/S1_R2.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix results/star/S1_
