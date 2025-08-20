# Step 3B-2: Gene-level counting with featureCounts

This step counts reads/fragments per gene from STAR-aligned BAM files.  
**Note:** Copy & paste into your terminal and replace placeholders `<>`.

## Command

```bash
#------------------------------------------------------------------------
featureCounts \
  -p \                         # paired-end mode (count fragments); remove for single-end
  -T <threads> \               # number of CPU threads
  -t exon \                    # feature type to count (exon)
  -g gene_id \                 # group features by 'gene_id' attribute in GTF
  -s <0|1|2> \                 # library strandedness: 0=unstranded, 1=stranded, 2=reverse
  -a <gene_annotation.gtf> \   # path to gene annotation (GTF)
  -o <gene_counts.txt> \       # output count table
  <aligned_bam_files>          # input BAM(s), e.g., results/star/*Aligned.sortedByCoord.out.bam
#------------------------------------------------------------------------
```

Notes

Paired-end data: keep -p. For single-end, remove -p.

Make sure the GTF matches the genome/annotation version used by STAR.

If you used name-sorted BAMs, add -B/-C options as needed (see featureCounts docs).

## Example
```
# Count fragments for multiple samples (paired-end, reverse-stranded libraries, 8 threads)
featureCounts \
  -p \
  -T 8 \
  -t exon \
  -g gene_id \
  -s 2 \
  -a ref/gencode.v46.annotation.gtf \
  -o results/counts/gene_counts.txt \
  results/star/*_Aligned.sortedByCoord.out.bam
```
