# Step7_ssGSEA_MSIGDB_from_counts.md

## Pipeline
`gene_counts.txt` → **DESeq2** VST → save normalized matrix → **GSVA::ssGSEA** with **MSigDB (e.g., HALLMARK or C7)** → export scores and heatmap of significant sets

## Inputs
- **gene_counts.txt**  
  Raw count matrix (rows = genes, cols = samples).  
- **sample.tsv**  
  Sample metadata file (must include sample identifiers).  
- **MSigDB gene sets**  
  Retrieved via `msigdbr` R package (choose collection = "H" for HALLMARK, "C7" for Immunologic Signatures, etc.).

## Outputs
- `07_input_normalized_vst.csv`  
  Normalized VST expression matrix.  
- `07-ssGSEA_scores.csv`  
  ssGSEA enrichment scores for each sample × gene set.  
- `07-ssGSEA_heatmap.png`  
  Heatmap of significant gene sets across samples.  

## Notes
- This script is general and can be applied to **any MSigDB collection**.  
- In this example, **HALLMARK** gene sets are used for demonstration, but you may substitute **C7** or other collections depending on your study.  
- Gene identifiers must match (HGNC symbols are typically required).  
