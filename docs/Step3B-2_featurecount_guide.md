# Step 3B-2. Gene-level Counting with featureCounts

## Tools 
- **featureCounts** (e.g., v2.0.3+)

## Inputs
- Aligned BAM files (from STAR)  
- Gene annotation GTF file (e.g., GENCODE, Ensembl)  

## Outputs
- `gene_counts.txt` : Gene-level count matrix (rows = genes, columns = samples)

## Key Options
- `-p` : Enable paired-end mode (count fragments instead of reads).  
- `-T <threads>` : Number of CPU threads.  
- `-t exon` : Count only exon features.  
- `-g gene_id` : Aggregate counts by the `gene_id` attribute in the GTF.  
- `-s <0|1|2>` : Library strandedness (0 = unstranded, 1 = stranded, 2 = reverse-stranded).  
- `-a <gene_annotation.gtf>` : Path to gene annotation file used for alignment.  
- `-o <gene_counts.txt>` : Output file containing the count matrix.  

## Notes
- The GTF annotation must match the genome reference used for STAR alignment.  
- For single-end libraries, omit `-p`.  
- For downstream DESeq2/edgeR analysis, `gene_counts.txt` serves as the raw input count matrix.
