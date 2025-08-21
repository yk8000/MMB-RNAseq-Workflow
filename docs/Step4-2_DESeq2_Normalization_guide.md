# Normalization with DESeq2

This step performs normalization of raw count data using **DESeq2**, producing variance-stabilized (VST) and regularized log-transformed (rlog) expression matrices for downstream analyses (e.g., PCA, clustering, differential expression).


## Option A. Salmon (tximport counts)

**Inputs**
- `samples.tsv` : sample metadata file  
- `gene_counts_salmon.csv` : gene-level count matrix summarized from Salmon (via tximport)  

**Outputs**
- `normalized_vst.csv` : variance-stabilized expression matrix  
- `normalized_rlog.csv` : rlog-transformed expression matrix  

**Notes**
- The design formula can be extended to include additional covariates, e.g. `~ batch + condition`.


## Option B. STAR + featureCounts

**Inputs**
- `samples.tsv` : sample metadata file  
- `gene_counts.txt` : output from featureCounts  

**Outputs**
- `normalized_vst.csv` : variance-stabilized expression matrix  
- `normalized_rlog.csv` : rlog-transformed expression matrix  

**Notes**
- Only the **Geneid** and raw count columns (`S1`, `S2`, …) are used.  
- Annotation columns (`Chr`, `Start`, `End`, `Strand`, `Length`) are ignored.  
- `sample_id` in `samples.tsv` must exactly match the column names in `gene_counts.txt`.



S3           Normal      Batch2
S4           Normal      Batch2
