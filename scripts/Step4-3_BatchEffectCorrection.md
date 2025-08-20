# Step 4-3. Batch Effect Correction (Optional)

Batch effects (e.g., sequencing run, library prep batch) can be addressed **during normalization** (by modeling `batch` in the DESeq2 design) or **after normalization** (by regressing out `batch` from a VST matrix). This guide provides both approaches in one place.

---

## Tools
- **R** (≥4.0)
- **DESeq2**
- **readr**, **dplyr**, **tibble**
- **limma** (for Option B)

---

## Inputs (common)
- `samples.tsv` (must contain at least: `sample_id`, `condition`, `batch`)
- **Option A**: a **raw count matrix** (either featureCounts `gene_counts.txt` or Salmon-derived `gene_counts_salmon.csv`)
- **Option B**: a **VST matrix** from Step 4-2 (e.g., `normalized_vst.csv`)

> Tip: `sample_id` values must match the count/vst matrix column names exactly.

---

## Option A — Model batch in the DESeq2 design (recommended when batch is known)

**When to use**
- You have raw counts and want DESeq2 to **jointly model** batch and condition in normalization and dispersion estimation.

**What you’ll get**
- VST and/or rlog matrices where batch effects are accounted for via the design formula.

```
# Rscript
# Option A: Account for batch in the DESeq2 design
# Design formula: ~ batch + condition

  library(readr)
  library(dplyr)
  library(tibble)
  library(DESeq2)

set.seed(1234)

# ---- Load sample metadata (must include: sample_id, condition, batch) ----
s <- read_tsv("samples.tsv", show_col_types = FALSE)

# ---- Load raw counts (choose ONE loader) ----

## (A1) featureCounts output (gene_counts.txt): keep Geneid + sample columns
fc <- read_tsv("gene_counts.txt", comment = "#", show_col_types = FALSE)
cm <- fc %>%
  select(Geneid, all_of(s$sample_id)) %>%   # keep only gene ID + sample columns
  column_to_rownames("Geneid") %>%          # gene IDs as row names
  as.matrix()

## (A2) Salmon-derived counts (gene_counts_salmon.csv, 1st column = gene IDs)
# cm <- read_csv("gene_counts_salmon.csv") %>%
#   column_to_rownames(1) %>%
#   as.matrix()

# ---- Align columns to samples.tsv and coerce to integer counts for DESeq2 ----
cm <- round(cm[, s$sample_id, drop = FALSE])

# ---- Build colData including batch, and create DESeq2 object ----
coldata <- data.frame(
  condition = s$condition,
  batch     = s$batch,
  row.names = s$sample_id
)

# Design models batch first, then condition of interest
dds <- DESeqDataSetFromMatrix(countData = cm, colData = coldata, design = ~ batch + condition)

# ---- Standard size-factor & dispersion estimation ----
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# ---- Variance-stabilized and rlog-transformed matrices (design-aware with blind = FALSE) ----
vst_mat  <- assay(vst(dds, blind = FALSE))
rlog_mat <- assay(rlog(dds, blind = FALSE))

# ---- Save results ----
write.csv(vst_mat,  "vst_design_batch_condition.csv")   # batch-adjusted via design
write.csv(rlog_mat, "rlog_design_batch_condition.csv")

message("[Batch Option A] Saved: vst_design_batch_condition.csv, rlog_design_batch_condition.csv")
```

