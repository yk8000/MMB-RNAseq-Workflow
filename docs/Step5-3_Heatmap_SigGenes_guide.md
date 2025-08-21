# Step5-3. Heatmap of Significant Genes

## Tools
- **DESeq2** (VST)
- **pheatmap**
- **org.Hs.eg.db**, **AnnotationDbi**
- **data.table**, **readr**

## Inputs
- `DEG_results_raw.csv` : DESeq2 results from **Step5-1** (includes `log2FoldChange`, `padj`)  
- `gene_counts.txt` : Raw gene-level count matrix (featureCounts or equivalent)  
- `samples.tsv` : Sample metadata (`condition` required; optional `patient_id`, `batch`)  

## Outputs
- `heatmap_sig_genes_2.png` : Clustered heatmap of significant genes across samples

## Key Options
- **Significance**: `padj < 0.05` and `|log2FC| > 1`  
- **Top genes cap**: `max_genes = 50` by smallest adjusted p-value  
- **Transformation**: VST (`varianceStabilizingTransformation`, `blind = TRUE`)  
- **Scaling**: Row-wise Z-score (per gene)  
- **Annotations**: `condition`, `patient_id`, `batch` shown on columns  
- **Clustering**: `cluster_rows = TRUE`, `cluster_cols = TRUE`  
- **Palette**: `colorRampPalette(c("navy","white","firebrick3"))(101)`  
- **Fonts**: `fontsize_row = 12`, `fontsize_col = 14`, `fontsize = 14`

## Notes
- Visualizes significant DEGs selected from **Step5-1** results; thresholds and the 50-gene cap are adjustable.  
- Samples are aligned between `gene_counts.txt` and `samples.tsv`; `condition` is treated as a factor with `"pre"` as the reference level.  
- Gene labels are mapped from Ensembl IDs to symbols via `org.Hs.eg.db` when available; otherwise Ensembl IDs are used.  
- If `batch` is absent in `samples.tsv`, a default single-level factor is assigned to avoid model errors.  
- The heatmap uses VST-transformed counts followed by **row-wise Z-score**, emphasizing relative expression patterns per gene across samples.  
- `set.seed(1234)` is used for reproducibility of clustering.  
