# Differential Expression Analysis with DESeq2
This step performs differential gene expression analysis using DESeq2 on raw count data. It compares biological conditions (e.g., pre vs post) and produces a results table containing log2 fold changes, p-values, and adjusted p-values. Example filtering is shown for padj < 0.05 & |log2FC| > 1.

## Inputs
samples.tsv : sample metadata file (must include condition, optional batch, patient_id)
gene_counts.txt : raw gene-level count matrix (featureCounts or Salmon summarized counts)

## Outputs
DEG_results_raw.csv : DESeq2 results for all genes
DEG_results_sig_padj0.05_log2FC1.csv : filtered DEG list (padj < 0.05 & |log2FC| > 1)

## Notes
Always use raw counts (not normalized values) as DESeq2 input.
The first level of condition (e.g., "pre") is treated as the reference in fold change calculation.
Choose design carefully:
Paired design: ~ patient_id + condition
Unpaired design: ~ batch + condition or ~ condition
Thresholds for significance (e.g., padj, log2FC) should be adjusted to study context.
Results will be used for downstream visualization and pathway analysis.
