# Step5-1. Differential Gene Expression (DGE) Analysis with DESeq2  
This step performs differential gene expression (DGE) analysis on raw count data using **DESeq2**, supporting both paired and unpaired designs, with optional batch adjustment.  

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
# - condition: biological group of interest (e.g., pre vs post treatment).
#              The first level ("pre") is treated as reference; log2FC is relative to this.
# - patient_id: subject identifier, used to model paired designs (pre/post from same patient).
# - batch: technical variable, used to account for processing or sequencing batch effects.
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

## Code Explanation
gene_counts.txt
```
gene_id         S1   S2   S3   S4
ENSG00000141510 120  98   150  130
ENSG00000153563 35   22   60   55
ENSG00000111640 5500 5100 6000 5800
ENSG00000075624 4800 4700 5200 5100
```

sample.tsv
```
sample_id   condition   patient_id   batch
P1-S1       pre         P1           B1
P1-S2       post        P1           B1
P2-S1       pre         P2           B2
P2-S2       post        P2           B2

```

- condition: Defines the biological groups being compared (e.g., pre vs post). The first level ("pre") is the reference; log2FoldChange values are computed relative to this reference.
-patient_id: Used in paired designs to explicitly model pre/post measurements from the same patient, thereby controlling for inter-patient variability.
- batch: Represents technical effects (e.g., sequencing run, library prep). Including batch in the design accounts for such unwanted variation.
- DESeqDataSetFromMatrix(...): Construct the DESeq2 dataset using the chosen design formula.
- DESeq(dds): Performs normalization, dispersion estimation, and hypothesis testing.
- results(dds, contrast = c("condition","post","pre")): Computes log2(post/pre). Positive values = higher in post; negative values = higher in pre.
- Example filtering with padj < 0.05 & |log2FC| > 1 is shown.
