# Step 7. Signature Analysis with ssGSEA
Using DESeq2 VST-normalized data (or TPM), we perform **single-sample GSEA (ssGSEA)** with MSigDB or custom gene sets to compute signature scores per sample and visualize them.  

## Pipeline
`gene_counts.txt` → **DESeq2 VST** *or* **TPM** → save normalized matrix → **GSVA::ssGSEA** → export results (score tables, heatmap)  

## Input data choice (important)
- **If analyzing only your own dataset** → use **VST** (DESeq2 variance-stabilized counts), recommended  
- **If integrating with external datasets (e.g., TCGA)** → use **TPM** for consistent scaling  

## Gene sets
- This script demonstrates **MSigDB HALLMARK** as an example  
- You may replace it with other MSigDB collections (e.g., **C7 Immunologic Signatures**) or your **custom gene sets**  
- For custom sets, the **GMT format** is recommended, which can be imported via `GSVA::getGmt("your_geneset.gmt")`  

## Inputs
- `gene_counts.txt`: gene count matrix (rows = genes, columns = samples); alternatively a TPM matrix can be used  
- `samples.tsv`: sample metadata (rows = sample IDs). Should contain at least `condition` (e.g., pre/post); `patient_id` and `batch` are optional but useful for annotation/statistics  

## Outputs
- `07_input_normalized_vst.csv`: normalized expression matrix (VST or TPM)  
- `07-ssGSEA_*_scores.csv`: ssGSEA score table (rows = gene sets, columns = samples)  
- `07-ssGSEA_sigsets_heatmap.png`: heatmap of significant signatures (row-wise Z-score, column annotation with condition, etc.)  

## Notes
- For Ensembl IDs, strip version suffix (e.g., `ENSG00000141510.12 → ENSG00000141510`) and map to HGNC SYMBOLs for compatibility with MSigDB  
- Statistical comparison is assumed for `pre` vs `post`: if sufficient pairs exist, a **paired Wilcoxon test** is applied; otherwise, an **unpaired Wilcoxon test** with **BH (FDR) correction** is used  
- Visualization focuses on significant signatures, displaying top sets with the strongest effect sizes for readability  
