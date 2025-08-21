# Heatmap of significant genes across samples
This script visualizes significant genes using DESeq2 results from Step5-1 (DEG_results_raw.csv)

```
# Rscript
library(DESeq2)
library(data.table)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)

set.seed(1234)  # for clustering stability

# ---------------------------
# 1) Inputs
# ---------------------------
# - res: DESeq2 results table (from Step5-1)
# - countData: raw counts (featureCounts or equivalent)
# - colData: sample metadata (must include condition; optional patient_id, batch)
res <- read.csv("DEG_results_raw.csv", row.names = 1, check.names = FALSE)
res <- as.data.frame(res)
res$padj[is.na(res$padj)] <- 1

countData <- read.table("gene_counts.txt", header = TRUE, row.names = 1, check.names = FALSE)
colData   <- read.table("samples.tsv",    header = TRUE, row.names = 1, check.names = FALSE, sep="\t")

# ---------------------------
# Align samples between counts and metadata
# ---------------------------
common_ids <- intersect(rownames(colData), colnames(countData))
countData  <- countData[, common_ids, drop = FALSE]
colData    <- colData[common_ids, , drop = FALSE]

# Define factors (reference level for condition = "pre")
colData$condition  <- factor(colData$condition, levels = c("pre","post"))
colData$patient_id <- factor(colData$patient_id)
if (!"batch" %in% colnames(colData)) colData$batch <- factor("B1") else colData$batch <- factor(colData$batch)

# ---------------------------
# 2) Select significant genes
# ---------------------------
# Thresholds: padj < 0.05 and |log2FC| > 1
alpha <- 0.05
lfc   <- 1
sig <- rownames(res)[which(res$padj < alpha & abs(res$log2FoldChange) > lfc)]

# Cap to top 50 by adjusted p-value (for readability)
max_genes <- 50
if (length(sig) > max_genes) {
  ord <- order(res[sig, "padj"], na.last = NA)
  sig <- sig[ord][1:max_genes]
}

# Keep only genes present in countData
sig <- intersect(sig, rownames(countData))
stopifnot(length(sig) > 1)

# ---------------------------
# 3) Transform counts (VST) and compute row-wise Z-scores
# ---------------------------
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countData),
                              colData   = colData,
                              design    = ~ condition)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
mat <- assay(vsd)[sig, , drop = FALSE]

# Row-wise Z-score
mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0

# ---------------------------
# 4) Map Ensembl IDs to gene symbols (if available)
# ---------------------------
ens_ids <- sub("\\..*$", "", rownames(mat_z)) # drop version suffix if any
syms <- mapIds(
  org.Hs.eg.db,
  keys     = ens_ids,
  keytype  = "ENSEMBL",
  column   = "SYMBOL",
  multiVals = "first"
)
# Replace rownames with SYMBOL if available
new_names <- ifelse(is.na(syms) | syms == "", rownames(mat_z), syms)
rownames(mat_z) <- new_names

# ---------------------------
# 5) Column annotations (sample-level)
# ---------------------------
ann_col <- data.frame(
  condition  = colData$condition,
  patient_id = colData$patient_id,
  batch      = colData$batch
)
rownames(ann_col) <- rownames(colData)

# ---------------------------
# 6) Plot & save heatmap
# ---------------------------
png("heatmap_sig_genes.png", width = 2400, height = 2400, res = 300)
pheatmap(
  mat_z,
  annotation_col = ann_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 8,
  fontsize = 14,
  scale = "none",
  color = colorRampPalette(c("navy","white","firebrick3"))(101),
  main = "Significant genes (vst, row Z-score)"
)
dev.off()

```

## Code explanation
gene_counts.txt: Gene-level raw count matrix (rows = genes, columns = samples).
```
gene_id         S1   S2   S3   S4
ENSG00000141510 120  98   150  130
ENSG00000153563 35   22   60   55
ENSG00000111640 5500 5100 6000 5800
ENSG00000075624 4800 4700 5200 5100
```
samples.tsv: Sample metadata; must include condition (and optional covariates such as batch, patient_id).
```
sample_id   condition   patient_id   batch
P1-S1          pre         P1           B1
P1-S2          post        P1           B1
P2-S1          pre         P2           B2
P2-S2          post        P2           B2
```
- res$padj[is.na(res$padj)] <- 1: Replace missing adjusted p-values with 1.
- intersect(...): Align sample IDs between counts and metadata tables.
- factor(...): Define condition, patient_id, and batch as factors (with "pre" as reference for condition).
- which(res$padj < alpha & abs(res$log2FoldChange) > lfc): Select significant genes (padj < 0.05, |log2FC| > 1).
- order(..., res[sig, "padj"]): Limit to top 50 significant genes by adjusted p-value.
- DESeqDataSetFromMatrix(...): Construct DESeq2 dataset for transformation.
- varianceStabilizingTransformation(...): Apply VST to normalize count data.
- scale(t(...)): Perform row-wise Z-score normalization.
- mapIds(..., org.Hs.eg.db): Map Ensembl IDs to gene symbols; replace row names when available.
- pheatmap(...): Draw clustered heatmap with row Z-scores, annotated by sample metadata, and save as PNG.
