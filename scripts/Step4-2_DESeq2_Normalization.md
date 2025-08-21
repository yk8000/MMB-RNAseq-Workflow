# Option A. Salmon (tximport counts)
## Normalization using DESeq2 on gene-level counts from Salmon (via tximport)

```
# Rscripts

library(readr)
library(tibble)
library(DESeq2)

set.seed(1234)

# ---------------------------
# 1) Load sample metadata
# ---------------------------
s <- read_tsv("samples.tsv", show_col_types = FALSE)

# ---------------------------
# 2) Load counts (Salmon-derived counts)
# ---------------------------
# gene_counts_salmon.csv: 1st column = gene IDs
cm <- read_csv("gene_counts_salmon.csv")
cm <- cm %>%
  column_to_rownames(1) %>%
  as.matrix()

# ---------------------------
# 3) Align columns with metadata and round to integers
# ---------------------------
cm <- round(cm[, s$sample_id, drop = FALSE])

# ---------------------------
# 4) Create colData (design: condition)
# ---------------------------
coldata <- data.frame(condition = s$condition,
                      row.names = s$sample_id)

# ---------------------------
# 5) Build DESeq2 object
# ---------------------------
dds <- DESeqDataSetFromMatrix(cm, coldata, ~ batch + condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# ---------------------------
# 6) Variance Stabilizing Transformation (VST)
# ---------------------------
vst <- vst(dds, blind = FALSE)
write.csv(assay(vst), "normalized_vst.csv")

# ---------------------------
# 7) Regularized log transformation (rlog)
# ---------------------------
rlg <- rlog(dds, blind = FALSE)
write.csv(assay(rlg), "normalized_rlog.csv")
")
```
# Option B. STAR + featureCounts
## Normalization using DESeq2 on gene-level counts from featureCounts
```
# Rscrips

library(readr)
library(dplyr)
library(tibble)
library(DESeq2)

# ---------------------------
# 1) Load sample metadata
# ---------------------------
s <- read_tsv("samples.tsv", show_col_types = FALSE)

# ---------------------------
# 2) Load counts (featureCounts output)
# ---------------------------
# Keep only Geneid + sample columns (drop annotation columns)
fc <- read_tsv("gene_counts.txt", comment = "#", show_col_types = FALSE)
cm <- fc %>%
  select(Geneid, all_of(s$sample_id)) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# ---------------------------
# 3) Align columns with metadata and round to integers
# ---------------------------
cm <- round(cm[, s$sample_id, drop = FALSE])

# ---------------------------
# 4) Create colData (design: condition)
# ---------------------------
coldata <- data.frame(condition = s$condition,
                      row.names = s$sample_id)

# ---------------------------
# 5) Build DESeq2 object and run normalization
# ---------------------------
dds <- DESeqDataSetFromMatrix(cm, coldata, ~ ~ batch + condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# ---------------------------
# 6) Variance Stabilizing Transformation (VST)
# ---------------------------
vst <- vst(dds, blind = FALSE)
write.csv(assay(vst), "normalized_vst.csv")

# ---------------------------
# 7) Regularized log transformation (rlog)
# ---------------------------
rlg <- rlog(dds, blind = FALSE)
write.csv(assay(rlg), "normalized_rlog.csv")

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
- sample_id: Sample name. Must match the column names in the count matrix.
- condition: Experimental group (e.g., Tumor, Normal). Used in the DESeq2 design formula.
- batch: Batch information (optional). If batch effect correction is required, specify the design formula as ~ batch + condition.

gene_counts.txt
```
Geneid	Chr	Start	End	Strand	Length	S1	S2	S3	S4
ENSG00000121410	chr1	11869	14409	+	1542	50	60	55	40
ENSG00000268895	chr1	14404	29570	-	16766	0	1	0	2
ENSG00000155657	chr1	34553	36081	+	1528	200	180	190	210
ENSG00000187634	chr1	69091	70008	-	917	10	12	8	9
```
- Geneid: Gene identifier (e.g., Ensembl).
- annotation columns: Chr, Start, End, Strand, Length provide annotation information.
- S1–S4: Raw counts for each sample. Must match the sample_id in samples.tsv.
