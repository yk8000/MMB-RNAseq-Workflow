# Step 5 Appendix. Mapping Ensembl Gene IDs to Gene Symbols

## Tools
- **R** (≥4.0)
- **org.Hs.eg.db**
- **AnnotationDbi**

---

## R Script

```r
#------------------------------------------------------------------------
library(org.Hs.eg.db)
library(AnnotationDbi)

# 1) Load DE results (rownames = Ensembl IDs)
res <- read.csv("DEG_results_raw.csv", row.names = 1, check.names = FALSE)

# 2) Drop version suffix from Ensembl IDs (e.g., ENSG00000141510.12 -> ENSG00000141510)
ens_ids <- sub("\\..*$", "", rownames(res))

# 3) Map Ensembl IDs to gene symbols
syms <- mapIds(
  org.Hs.eg.db,
  keys     = ens_ids,
  keytype  = "ENSEMBL",
  column   = "SYMBOL",
  multiVals = "first"
)

# 4) Add SYMBOL column
res$SYMBOL <- syms

# 5) Save
write.csv(res, "DEG_results_with_symbol.csv")
#------------------------------------------------------------------------
