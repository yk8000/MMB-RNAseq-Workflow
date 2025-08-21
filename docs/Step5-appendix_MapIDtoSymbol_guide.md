# Step5 Appendix. Mapping Ensembl Gene IDs to Gene Symbols
This appendix adds human gene symbols to DESeq2 results by mapping Ensembl gene IDs (optionally version-stripped) using Bioconductor annotations.

## Tools
- **R** (≥4.0)
- **org.Hs.eg.db**
- **AnnotationDbi**

## Inputs
- `DEG_results_raw.csv` : DESeq2 results from **Step5-1** (row names = Ensembl gene IDs)

## Outputs
- `DEG_results_with_symbol.csv` : DESeq2 results with an added `SYMBOL` column

## Key Options
- Strip version suffix from Ensembl IDs (e.g., `ENSG00000141510.12` → `ENSG00000141510`) before mapping.  
- Use `AnnotationDbi::mapIds(...)` with `keytype = "ENSEMBL"`, `column = "SYMBOL"`, `multiVals = "first"`.

## Notes
- Ensure the species-specific OrgDb matches your data (e.g., `org.Hs.eg.db` for human; `org.Mm.eg.db` for mouse).  
- If a symbol is not found, the value remains `NA`; consider keeping Ensembl IDs as fallback labels.  
- Keep genome/annotation releases consistent to avoid missing or outdated mappings.  
- Many-to-one mappings exist; `multiVals = "first"` takes the first match. If you need all mappings, consider alternatives (e.g., `multiVals = "CharacterList"`).
