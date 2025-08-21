# Volcano plot (EnhancedVolcano style)
This script visualizes the DESeq2 results from Step5-1 (DEG_results_raw.csv) 
by drawing a volcano plot highlighting significant up- and down-regulated genes.

```
# Rscript

library(EnhancedVolcano)

# ---------------------------
# 1) Load DESeq2 results (row names = genes/Ensembl IDs)
# ---------------------------
# - Input: DEG_results_raw.csv produced by DESeq2 (Step5-1)
# - row.names = gene identifiers (Ensembl IDs by default)
res <- read.csv("DEG_results_raw.csv", row.names = 1, check.names = FALSE)

# ---------------------------
# 2) Replace NA padj with 1 to keep genes in the plot
# ---------------------------
# - Some genes have NA for adjusted p-values (low counts etc.)
# - Replace with 1 so they appear as non-significant in the plot
res$padj[is.na(res$padj)] <- 1

# ---------------------------
# 3) Prepare labels (default = Ensembl IDs; map to gene symbols if available)
# ---------------------------
# - Default labels = Ensembl IDs
# - If org.Hs.eg.db & AnnotationDbi are available:
#     drop version suffix from Ensembl IDs (e.g., ENSG00000141510.15 → ENSG00000141510)
#     map to SYMBOL (HGNC gene symbol)
# - Fallback: keep Ensembl IDs if mapping fails
lab <- rownames(res)

if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)) {
  ens_ids <- sub("\\..*$", "", rownames(res)) # drop version suffix
  syms <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys     = ens_ids,
    keytype  = "ENSEMBL",
    column   = "SYMBOL",
    multiVals = "first"
  )
  lab <- ifelse(is.na(syms) | syms == "", rownames(res), syms)
}

# ---------------------------
# 4) Define significance thresholds and select top labels
# ---------------------------
# - alpha = FDR threshold (padj < 0.05)
# - lfc   = effect size threshold (|log2FC| > 1)
# - up_idx   = significant up-regulated genes
# - down_idx = significant down-regulated genes
# - Among each set, order by padj and select top 10 for labeling
alpha <- 0.05
lfc   <- 1

up_idx   <- which(res$padj < alpha & res$log2FoldChange >  lfc)
down_idx <- which(res$padj < alpha & res$log2FoldChange < -lfc)

up_ord   <- order(res$padj[up_idx], na.last = NA)
down_ord <- order(res$padj[down_idx], na.last = NA)

up_sel   <- lab[up_idx][up_ord][1:min(10, length(up_idx))]
down_sel <- lab[down_idx][down_ord][1:min(10, length(down_idx))]

# - selectLab = genes to highlight in plot
selectLab <- c(up_sel, down_sel)

# ---------------------------
# 5) Plot and save volcano plot
# ---------------------------
# - Output: volcano_pre_vs_post.png
# - EnhancedVolcano parameters:
#     x: log2FoldChange (effect size)
#     y: padj (adjusted p-value)
#     pCutoff: significance threshold (alpha)
#     FCcutoff: log2 fold change cutoff
#     selectLab: top genes to label
#     col: custom color scheme (grey = ns, green = low effect, blue = one side, red = sig)
#     labSize, axisLabSize, titleLabSize: adjust font sizes
#     drawConnectors: link labels to points
png("volcano_pre_vs_post.png", width = 3000, height = 1800, res = 300)
print(
  EnhancedVolcano(
    toptable    = res,
    lab         = lab,
    x           = "log2FoldChange",
    y           = "padj",
    pCutoff     = alpha,
    FCcutoff    = lfc,
    title       = "Volcano plot",
    # subtitle  = "EnhancedVolcano",
    xlab        = bquote(~Log[2]~"fold change (post / pre)"),
    ylab        = bquote(-Log[10]~"(adjusted p-value)"),
    caption     = paste("total =", nrow(res), "variables"),
    legendPosition = "right",
    col         = c("grey50","forestgreen","royalblue","red2"),
    pointSize   = 1.4,
    labSize     = 6,
    titleLabSize = 20,            
    subtitleLabSize = 16,         
    captionLabSize = 14,          
    axisLabSize = 16,             
    drawConnectors = TRUE,
    widthConnectors = 0.6,
    selectLab   = selectLab
  )
)
dev.off()

```
## Code Explanation
DEG_results_raw.csv (produced by Step5-1)
```
,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj
ENSG00000141510,102.5, 1.35,0.42, 3.20,0.0014,0.012
ENSG00000153563, 45.3,-2.10,0.55,-3.82,0.00013,0.004
ENSG00000111640,5300.7, 0.05,0.10, 0.50,0.61,0.89
ENSG00000075624,4900.2,-0.20,0.08,-2.50,0.012,0.20
ENSG00000123415,  75.6, 2.25,0.70, 3.21,0.0013,0.011
```

- baseMean: average normalized counts
- log2FoldChange: log2(post/pre) fold change
- padj: adjusted p-value (FDR) used for significance cutoff

