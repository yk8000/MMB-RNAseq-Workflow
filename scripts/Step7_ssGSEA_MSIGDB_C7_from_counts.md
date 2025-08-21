# Step 7. Immune Signature Analysis with ssGSEA  
This step computes **immune signature scores** using **ssGSEA (GSVA)** on a DESeq2 VST-normalized matrix built directly from raw counts.  
Pipeline: `gene_counts.txt` → **DESeq2** VST → save normalized matrix → **GSVA::ssGSEA** with **MSigDB Hallmark (H collection, example)** → export scores and a heatmap of significant sets.


```
# Rscript

library(DESeq2)
library(GSVA)
library(msigdbr)
library(pheatmap)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

set.seed(1234) # project convention

# ---------------------------
# 1) Load inputs
# ---------------------------
# gene_counts.txt: tab-delimited, header present, first column = gene IDs
counts <- read.table("gene_counts.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
counts <- as.matrix(counts)

# sample.tsv: tab-delimited, header present, first column = sample IDs (must match column names of counts)
colData <- read.table("samples.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Align samples between counts and metadata
common_samples <- intersect(colnames(counts), rownames(colData))
if (length(common_samples) == 0) stop("No shared sample IDs between counts (columns) and sample.tsv (rownames).")
counts  <- counts[, common_samples, drop = FALSE]
colData <- colData[common_samples, , drop = FALSE]

# If no 'condition' column exists, create a dummy intercept-only design
if (!"condition" %in% colnames(colData)) {
  message("No 'condition' column found in sample.tsv. Creating a dummy intercept-only design.")
  colData$condition <- factor("all")
  design_formula <- ~ 1
} else {
  colData$condition <- factor(colData$condition)
  design_formula <- ~ condition
}

# ---------------------------
# 2) Build DESeq2 object and VST
# ---------------------------
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = design_formula)
dds <- dds[rowSums(counts(dds)) > 0, ]  # keep genes with non-zero counts
vsd <- vst(dds, blind = TRUE)
expr_vst <- assay(vsd)  # genes x samples

# ---------------------------
# 3) Optional: convert Ensembl IDs to HGNC symbols (robust base R version)
#     - Strip version suffix (e.g., ENSG00000141510.12 -> ENSG00000141510)
#     - Map to SYMBOL; if unmapped, keep original ID
#     - Collapse duplicates by SYMBOL (average across rows)
# ---------------------------
looks_like_ensembl <- grepl("^ENSG\\d+", rownames(expr_vst))
if (any(looks_like_ensembl)) {
  ensembl_novers <- sub("\\..*$", "", rownames(expr_vst))
  syms <- mapIds(org.Hs.eg.db,
                 keys = ensembl_novers,
                 keytype = "ENSEMBL",
                 column = "SYMBOL",
                 multiVals = "first")
  syms_filled <- ifelse(is.na(syms) | syms == "", rownames(expr_vst), syms)

  # Check for duplicated sample names; make them unique to avoid downstream issues
  if (any(duplicated(colnames(expr_vst)))) {
    warning("Duplicated sample names detected; making them unique via make.unique().")
    colnames(expr_vst) <- make.unique(colnames(expr_vst))
  }

  # Build a data.frame and aggregate duplicates by SYMBOL using mean
  expr_df <- as.data.frame(expr_vst, check.names = FALSE)
  expr_df$SYMBOL <- syms_filled

  expr_agg <- aggregate(. ~ SYMBOL, data = expr_df, FUN = mean)

  rownames(expr_agg) <- expr_agg$SYMBOL
  expr_agg$SYMBOL <- NULL
  expr_vst <- as.matrix(expr_agg)
}

# ---------------------------
# 4) Save normalized matrix for Step 7
# ---------------------------
write.csv(expr_vst, "07_input_normalized_vst.csv", quote = FALSE)

# ---------------------------
# 5) Prepare Hallmark gene sets (MSigDB H: Hallmark, human, example)
# ---------------------------
h_all <- msigdbr(species = "Homo sapiens", category = "H")

# - Each row corresponds to one gene in one Hallmark set
# - Select columns for set name and gene symbol
hallmark <- h_all[, c("gs_name", "gene_symbol")]

# Convert to a named list (genes per gene set)
split_by_set <- split(hallmark$gene_symbol, hallmark$gs_name)
genesets <- lapply(split_by_set, unique)

# Restrict to genes present in the expression matrix
genesets <- lapply(genesets, function(g) intersect(g, rownames(expr_vst)))
genesets <- genesets[vapply(genesets, length, 1L) > 0]

if (length(genesets) == 0) stop("No Hallmark gene sets overlap with the expression matrix row names.")

# ---------------------------
# 6) ssGSEA (GSVA)
# ---------------------------
param  <- ssgseaParam(expr_vst, genesets, normalize = TRUE)
scores <- gsva(param)   

# ---------------------------
# 7) Save scores and heatmap
# ---------------------------
write.csv(scores, "07-ssGSEA_HM_scores.csv", quote = FALSE)

# Per-set Wilcoxon test: explicitly separate paired vs unpaired comparisons
test_sets_split <- function(scores, meta, min_pairs = 3) {
  stopifnot(all(colnames(scores) %in% rownames(meta)))
  meta <- meta[colnames(scores), , drop = FALSE]

  res <- lapply(rownames(scores), function(gs){
    df <- data.frame(
      score      = as.numeric(scores[gs, ]),  # ssGSEA score per sample
      patient_id = meta$patient_id,           # used to detect pairing
      condition  = meta$condition,            # pre/post grouping
      stringsAsFactors = FALSE
    )

    # Check if both pre/post samples exist for each patient -> if yes, do paired test
    pts_pre  <- unique(df$patient_id[df$condition == "pre"])
    pts_post <- unique(df$patient_id[df$condition == "post"])
    pts      <- intersect(pts_pre, pts_post)

    if (length(pts) >= min_pairs) {
      ## ----- paired: compute per-patient diffs (post - pre) and run signed-rank test -----
      diffs <- sapply(pts, function(id){
        pre_val  <- df$score[df$condition=="pre"  & df$patient_id==id][1]
        post_val <- df$score[df$condition=="post" & df$patient_id==id][1]
        post_val - pre_val
      })
      diffs <- diffs[is.finite(diffs)]
      stat  <- suppressWarnings(wilcox.test(diffs))  # Wilcoxon signed-rank test on diffs
      # Reference means (reported alongside)
      mean_pre  <- mean(df$score[df$condition=="pre"],  na.rm=TRUE)
      mean_post <- mean(df$score[df$condition=="post"], na.rm=TRUE)
      paired <- TRUE

    } else {
      ## ----- unpaired: fall back to rank-sum test between groups -----
      pre_vals  <- df$score[df$condition=="pre"]
      post_vals <- df$score[df$condition=="post"]
      stat  <- suppressWarnings(wilcox.test(post_vals, pre_vals, paired = FALSE))
      mean_pre  <- mean(pre_vals,  na.rm=TRUE)
      mean_post <- mean(post_vals, na.rm=TRUE)
      paired <- FALSE
    }

    data.frame(
      set       = gs,
      mean_pre  = mean_pre,
      mean_post = mean_post,
      diff      = mean_post - mean_pre,  # effect size reported
      p         = unname(stat$p.value),
      paired    = paired,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, res)
  out$FDR <- p.adjust(out$p, method = "BH")  # multiple testing correction
  out
}

# Run tests and filter significant sets
res <- test_sets_split(scores, colData)
sig <- subset(res, FDR < 0.05)                # FDR threshold
sig <- subset(sig, abs(diff) >= 0.05)         # optional effect-size filter

# If too many, keep top-N by |diff|
N <- 50
sig <- sig[order(abs(sig$diff), decreasing = TRUE), , drop = FALSE]
sig_top <- head(sig, N)

# Subset matrix for heatmap
if (nrow(sig_top) == 0) stop("No significant gene sets at the chosen thresholds.")
mat <- scores[sig_top$set, , drop = FALSE]

# Column annotations (optional, shown if columns exist in metadata)
ann <- data.frame(
  batch      = colData[colnames(mat), "batch"],
  patient_id = colData[colnames(mat), "patient_id"],
  condition  = colData[colnames(mat), "condition"]
)
rownames(ann) <- colnames(mat)

# Plot heatmap (row-wise Z-score via scale(); annotations show batch/patient/condition)
png("07-ssGSEA_sigsets_heatmap.png", width = 2800, height = 2000, res = 300)
pheatmap(
  scale(mat),
  annotation_col = ann,
  color  = colorRampPalette(c("navy","white","firebrick3"))(100),
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 45,
  border_color = NA,
  fontsize_row = 6,
  fontsize_col = 9,
  main = "Significant immunologic signatures (Wilcoxon, FDR<0.05)"
)
dev.off()
```

## Code Explanation

- **alpha <- 0.05; lfc <- 1**  
  Define thresholds for significance (adjusted p-value < 0.05 and |log2FC| > 1).
- **res$padj[is.na(res$padj)] <- 1**  
  Replace missing adjusted p-values with 1, so they are treated as non-significant.
- **up_ids_raw / down_ids_raw**  
  Split DEGs into UP (log2FC > 1) and DOWN (log2FC < -1) sets under the padj threshold.
- **strip_version <- function(x) sub("\\..*$", "", x)**  
  Remove Ensembl version suffix (e.g., ENSG00000141510.12 → ENSG00000141510).
- **bitr(...)**  
  Map Ensembl IDs to ENTREZ IDs and SYMBOLs for KEGG analysis. Drop unmapped IDs.
- **run_enrich_go_ens(...)**  
  Perform GO enrichment using Ensembl IDs as input; results are converted to readable SYMBOLs.
- **simplify_bp(...)**  
  Reduce redundancy in GO Biological Process terms by clustering similar terms.
- **run_enrich_kegg(...)**  
  Perform KEGG pathway enrichment using ENTREZ IDs.
- **one_side(...)**  
  Wrapper to handle one DEG side (UP or DOWN):  
  1. Save DEG lists (Ensembl + mapped ENTREZ)  
  2. Run GO (BP/CC/MF) and KEGG enrichment  
  3. Save results as CSVs and dot plots  
  4. Handle empty sets gracefully by outputting empty CSVs.
