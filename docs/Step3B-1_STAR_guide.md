# Step 3B-1. Genome Indexing and Read Alignment with STAR

**Purpose**  
Build a genome index and align RNA-seq reads to the reference genome using STAR.  
This step produces coordinate-sorted BAM files for downstream quantification.

**Tools**  
- STAR (e.g., v2.7+)

**Inputs**  
- Reference genome FASTA (e.g., `GRCh38.primary_assembly.genome.fa`)  
- Gene annotation GTF (e.g., `gencode.v46.annotation.gtf`)  
- Raw FASTQ files (paired-end or single-end)  

**Outputs**  
- `ref/STAR_index/` : STAR genome index directory  
- `results/star/<sample>_Aligned.sortedByCoord.out.bam` : Sorted BAM file per sample  
- `results/star/<sample>_Log.final.out` : Alignment summary statistics  
- `results/star/<sample>_SJ.out.tab` : Detected splice junctions  

**Notes**  
- For gzipped FASTQ files (`*.fastq.gz`), add the option: `--readFilesCommand zcat`.  
- Use the same reference genome and annotation version consistently across all steps.  
- Ensure sufficient disk space: the STAR index can be >30 GB for human genome.
