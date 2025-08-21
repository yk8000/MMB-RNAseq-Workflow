# Step8. Immune cell deconvolution and cell state inference

Estimate immune/stromal cell composition and infer cell states/ecotypes from bulk RNA‑seq using web-based tools (no script needed for this step).


## Input data
- **CIBERSORTx, xCell**: TPM matrix  
- **MCP-counter**: Raw gene-level count matrix  
- **Ecotyper**: TPM or normalized count matrix (recommended: TPM)

## Tools
- **CIBERSORTx** (https://cibersortx.stanford.edu/)  
  → Estimates relative immune cell fractions from bulk expression profiles.  
- **MCP-counter** (https://github.com/ebecht/MCPcounter)  
  → Provides abundance scores for immune and stromal populations.  
- **xCell** (https://xcell.ucsf.edu/)  
  → Performs enrichment analysis for 64 immune and stromal cell types.  
- **Ecotyper** (https://ecotyper.stanford.edu/)  
  → Infers cell states and ecotypes from gene expression data.


## Output
- Immune cell proportions for each sample  
- Cell states and ecotypes inferred by **Ecotyper**


## Notes
- Each tool has specific input requirements (see above).  
- TPM is generally recommended when using **multiple tools** for consistency.  
- Results can be integrated with DEGs and ssGSEA outputs (Step6–7) for downstream analyses.
