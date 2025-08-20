## A)	DESeq2 design: include batch factor (from raw counts)
Normalization including batch factor in the design (~ batch + condition)

```
# Rscripts

library(readr)
library(dplyr)
library(tibble)
library(DESeq2)

# 1) Load sample metadata (samples.tsv: must contain sample_id, condition, batch)
s <- read_tsv("samples.tsv", show_col_types = FALSE)

# 2) Load counts (choose one loader)
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

# 3) Align columns with metadata and round to integers
cm <- round(cm[, s$sample_id, drop = FALSE])

# 4) Create colData including batch
coldata <- data.frame(condition = s$condition,
                      batch     = s$batch,
                      row.names = s$sample_id)

# 5) Create DESeq2 dataset with batch factor in design (~ batch + condition)
dds <- DESeqDataSetFromMatrix(cm, coldata, ~ batch + condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# 6) Variance Stabilizing Transformation (VST), design-aware
vst_mat <- assay(vst(dds, blind = FALSE))

# 7) Save VST matrix
write.csv(vst_mat, "vst_design_batch_condition.csv")
```
## Option B. limma::removeBatchEffect
Batch effect removal from an already normalized VST matrix
