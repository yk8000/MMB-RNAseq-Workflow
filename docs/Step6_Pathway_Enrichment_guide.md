# Step6. Pathway and Functional Enrichment Analysis with clusterProfiler
This step performs **Over-Representation Analysis (ORA)** using `clusterProfiler`.  
From the DESeq2 results (`DEG_results_raw.csv`) generated in Step5-1, the script:  
- Defines UP- and DOWN-regulated gene sets (based on adjusted p-value and log2FC thresholds).  
- Runs **GO enrichment** (Biological Process, Cellular Component, Molecular Function) using ENSEMBL IDs.  
- Runs **KEGG enrichment** using ENTREZ IDs mapped from Ensembl.  
- Outputs results as **CSV tables** and **dot plots (PNG, 300 dpi)** into `./UP` and `./DOWN` directories.  


## Inputs
- **DEG_results_raw.csv**  
  - A CSV file from Step5-1 with DESeq2 results.  
  - Row names = Ensembl gene IDs (may include version suffix).  
  - Must contain at least:  
    - `log2FoldChange`  
    - `padj` (adjusted p-value)


## Outputs
- `UP/` and `DOWN/` directories containing:  
  - `GO_BP_*.csv`, `GO_CC_*.csv`, `GO_MF_*.csv`: tables of enriched GO terms  
  - `KEGG_*.csv`: tables of enriched KEGG pathways  
  - Corresponding dot plots: `*_dotplot.png`  
  - DEG lists used for enrichment:  
    - `*_genes_used_GO_ENSEMBL.csv` (Ensembl IDs)  
    - `*_genes_mapped_ENTREZ.csv` (mapped ENTREZ IDs)  


## Notes
- **GO analysis**:  
  - Uses `keyType="ENSEMBL"`.  
  - `readable=TRUE` converts IDs to SYMBOL in the results.  
  - Universe = all tested Ensembl IDs.  
- **KEGG analysis**:  
  - Requires ENTREZ IDs, obtained via `clusterProfiler::bitr`.  
  - Universe = all successfully mapped ENTREZ IDs.  
- **Parameters**:  
  - `alpha = 0.05` (FDR threshold for DEGs)  
  - `lfc = 1` (absolute log2FC cutoff = 2-fold change)  
  - Dot plots show top 10 categories by default (`topn_go_bp`, `topn_go_etc`, `topn_kegg`).  
- **Reproducibility**:  
  - `set.seed(1234)` is included for consistency across project scripts (though enrichment itself is deterministic).  

