# Step 5 Appendix. Mapping Ensembl Gene IDs to Gene Symbols
This step adds gene symbols to the DESeq2 results by mapping Ensembl gene IDs (optionally version-stripped) to human gene symbols using the `org.Hs.eg.db` annotation database.


```
# R Script

# ---------------------------
# 1) Inputs
# ---------------------------
# - res: DESeq2 results table (from Step5-1)
res <- read.csv("DEG_results_raw.csv", row.names = 1, check.names = FALSE)

# ---------------------------
# 2) Drop version suffix from Ensembl IDs
#    (e.g., ENSG00000141510.12 -> ENSG00000141510)
# ---------------------------
ens_ids <- sub("\\..*$", "", rownames(res))

# ---------------------------
# 3) Map Ensembl IDs to gene symbols
# ---------------------------
syms <- mapIds(
  org.Hs.eg.db,
  keys     = ens_ids,
  keytype  = "ENSEMBL",
  column   = "SYMBOL",
  multiVals = "first"
)

# ---------------------------
# 4) Add SYMBOL column
# ---------------------------
res$SYMBOL <- syms

# ---------------------------
# 5) Save
# ---------------------------
write.csv(res, "DEG_results_with_symbol.csv")
