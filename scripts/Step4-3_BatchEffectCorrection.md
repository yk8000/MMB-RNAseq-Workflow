# Step 4-3: Batch Effect Correction
## A)	DESeq2 design: include batch factor (from raw counts)
In this approach, batch effects are addressed directly within the **DESeq2** model design. By specifying a formula such as `~ batch + condition`, DESeq2 accounts for both biological condition and batch simultaneously during normalization and differential expression analysis. This method ensures that unwanted technical variation is adjusted while retaining true biological signals.  

```
# Rscript

library(readr)
library(dplyr)
library(tibble)
library(DESeq2)

# ---------------------------
# 1) Load sample metadata
# ---------------------------
# samples.tsv: must contain sample_id, condition, batch
s <- read_tsv("samples.tsv", show_col_types = FALSE)

# ---------------------------
# 2) Load counts (choose one loader)
# ---------------------------
## 2-1) featureCounts output (gene_counts.txt)
fc <- read_tsv("gene_counts.txt", comment = "#", show_col_types = FALSE)
cm <- fc %>%
  select(Geneid, all_of(s$sample_id)) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

## 2-2) (Alternative) Salmon-derived counts (gene_counts_salmon.csv)
# cm <- read_csv("gene_counts_salmon.csv") %>%
#          column_to_rownames(1) %>%
#          as.matrix()

# ---------------------------
# 3) Align columns with metadata and round to integers
# ---------------------------
cm <- round(cm[, s$sample_id, drop = FALSE])

# ---------------------------
# 4) Create colData including batch
# ---------------------------
coldata <- data.frame(condition = s$condition,
                      batch     = s$batch,
                      row.names = s$sample_id)

# ---------------------------
# 5) Create DESeq2 dataset with batch factor in design (~ batch + condition)
# ---------------------------
# - 'batch' : specify batch information (from samples.tsv)
# - 'design': model.matrix(~ condition) keeps condition effects intact
dds <- DESeqDataSetFromMatrix(cm, coldata, ~ batch + condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# ---------------------------
# 6) Variance Stabilizing Transformation (VST), design-aware
# ---------------------------
vst_mat <- assay(vst(dds, blind = FALSE))

# ---------------------------
# 7) Save VST matrix
# ---------------------------
write.csv(vst_mat, "vst_design_batch_condition.csv")

```
## Option B. limma::removeBatchEffect
Batch effect removal from an already normalized VST matrix

```
# Rscript

library(readr)
library(limma)

# ---------------------------
# 1) Load sample metadata
# ---------------------------
s <- read_tsv("samples.tsv", show_col_types = FALSE)

# ---------------------------
# 2) Load VST matrix from previous step
# ---------------------------
# - normalized_vst.csv produced by DESeq2 (rows = genes, cols = samples)
vst_mat <- as.matrix(read.csv("normalized_vst.csv", 
                              row.names = 1, 
                              check.names = FALSE))

# ---------------------------
# 3) Align columns with metadata
# ---------------------------
# - ensure order of columns in vst_mat matches s$sample_id
vst_mat <- vst_mat[, s$sample_id, drop = FALSE]

# ---------------------------
# 4) Remove batch effect using limma::removeBatchEffect
# ---------------------------
# - 'batch' : specify batch information (from samples.tsv)
# - 'design': model.matrix(~ condition) keeps condition effects intact
design <- model.matrix(~ s$condition)
vst_bc <- removeBatchEffect(vst_mat, batch = s$batch, design = design)

# ---------------------------
# 5) Save batch-corrected VST matrix
# ---------------------------
# - output file: vst_batchCorrected_limma.csv
write.csv(vst_bc, "vst_batchCorrected_limma.csv")

```

## Code Explanation
samples.tsv
```
sample_id	condition	batch
S1	Tumor	Batch1
S2	Tumor	Batch1
S3	Normal	Batch2
S4	Normal	Batch2
```

## Notes
- Use **Option A (design-based)** if you are still at the raw-count stage and plan to run differential expression (`~ batch + condition`).  
- Use **Option B (limma post-hoc)** when you already have a normalized VST matrix and want to remove batch effects for visualization (PCA, heatmaps, clustering).  
