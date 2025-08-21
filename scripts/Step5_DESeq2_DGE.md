# Differential expression on raw counts using DESeq2

```
# Rscript

library(DESeq2)

set.seed(1234)  # for reproducibility

# ---------------------------
# 1) Load inputs
# ---------------------------
# - gene_counts.txt: rows = genes, cols = samples (raw counts)
# - samples.tsv    : rows = samples; must include 'condition'
countData <- read.table("gene_counts.txt", header = TRUE, row.names = 1, check.names = FALSE)
colData   <- read.table("samples.tsv",   header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")

# ---------------------------
# 2) Align samples
# ---------------------------
common_ids <- intersect(rownames(colData), colnames(countData))
countData  <- countData[, common_ids, drop = FALSE]
colData    <- colData[common_ids, , drop = FALSE]

# ---------------------------
# 3) Factors
# ---------------------------
colData$condition  <- factor(colData$condition, levels = c("pre","post"))
colData$patient_id <- factor(colData$patient_id)
colData$batch      <- factor(colData$batch)

# ---------------------------
# 4) Build DESeq2 object
#    Choose ONE of the following designs.
# ---------------------------

## 4-1) Paired design (same patients pre/post)
##     If batches also exist across the same patients, use: ~ batch + patient_id + condition
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(countData),
  colData   = colData,
  design    = ~ patient_id + condition
)

## 4-2) Unpaired design
##     If batch exists: design = ~ batch + condition
##     If no batch    : design = ~ condition
# dds <- DESeqDataSetFromMatrix(
#   countData = as.matrix(countData),
#   colData   = colData,
#   design    = ~ batch + condition
# )

# ---------------------------
# 5) Run DESeq2
# ---------------------------
dds <- DESeq(dds)

# ---------------------------
# 6) Contrast: post vs pre
#     - contrast = c("condition", "post", "pre") → log2(post/pre)
#     - log2FoldChange > 0 → higher in post
#     - log2FoldChange < 0 → higher in pre
# ---------------------------
res <- results(dds, contrast = c("condition", "post", "pre"))

# ---------------------------
# 7) Save raw results
# ---------------------------
write.csv(as.data.frame(res), "DEG_results_raw.csv", row.names = TRUE)

# ---------------------------
# 8) Example filtering
#     padj < 0.05 & |log2FC| > 1
# ---------------------------
res_df <- as.data.frame(res)
res_df$padj[is.na(res_df$padj)] <- 1
sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig, "DEG_results_sig_padj0.05_log2FC1.csv", row.names = TRUE)
```
