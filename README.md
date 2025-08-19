# Bulk RNA-seq Workflow for Immuno-Oncology

This repository provides a **bulk RNA-seq analysis pipeline** applicable to both public and in-house datasets.  
Target audience: cancer immunology and clinical researchers.

## Scope
- Covers preprocessing, normalization, differential expression analysis, pathway enrichment (GO/KEGG), immune signature analysis (ssGSEA/GSVA), and immune cell deconvolution (CIBERSORTx, etc.).
- Provides scripts and execution order to regenerate figures and tables (**no raw data included**).

## Repository Structure
```
project/
├─ README.md
├─ LICENSE # MIT (code)
├─ scripts/
│  ├─ Step1-1_<desc>.R
│  ├─ Step1-2_<desc>.R
│  ├─ Step2-1_<desc>.R
│  └─ ...
└─ docs/ # Supplementary notes
```

## Naming & Conventions
- **Scripts**: `StepX-Y_<short-description>.R` (e.g., `Step2-1_ssGSEA_from_counts.R`)
- **Outputs**: start with the corresponding step number (e.g., `Step2-1_ssGSEA_scores.csv`)
- Figures are saved in `png()` format (resolution = 300 dpi).

## Data Policy
- This repository does **not** contain data (`data/` is not provided).
- External resources (e.g., MSigDB gene sets, CIBERSORT signatures) are **not redistributed**; links are provided instead.
- Restricted-access datasets (e.g., EGA) are not included. Only instructions for access are provided.

### External Resources (Examples)
- MSigDB: <https://www.gsea-msigdb.org/gsea/msigdb/>
- CIBERSORTx: <https://cibersortx.stanford.edu/>
- xCell: <https://aran-lab.com/xcell2-vignette/>

## How to Run (Overview)

### Step 1–2: Quality Control and Preprocessing
1. **Quality Check**  
   - Tools: FastQC, MultiQC  
   - Input: Raw FASTQ files  
   - Output: Quality reports  
2. **Adapter Trimming & Filtering**  
   - Tools: Trim Galore! (uses Cutadapt internally)  
   - Input: Raw FASTQ files  
   - Output: Cleaned FASTQ files  

### Step 3–4: Expression Matrix Generation
3. **Alignment & Quantification**  
   - Option A (alignment-free): Salmon → tximport (R)  
     - Input: Cleaned FASTQ, Reference transcriptome (FASTA)  
     - Output: Raw gene-level count matrix  
   - Option B (alignment-based): STAR → featureCounts  
     - Input: Cleaned FASTQ, Reference genome (FASTA), Annotation (GTF)  
     - Output: Raw gene-level count matrix  
4. **Normalization & Batch Correction**  
   - TPM calculation → for deconvolution  
   - DESeq2 vst/rlog → for GSVA and visualization  
   - Optional: batch correction with limma or DESeq2 design  
   - Outputs: TPM matrix, Normalized count matrix, Batch-corrected matrix  

### Step 5–6: Differential Expression and Pathway Analysis
5. **Differential Gene Expression Analysis (DGEA)**  
   - Tool: DESeq2  
   - Input: Raw count matrix + sample metadata  
   - Output: Differentially Expressed Genes (DEGs)  
6. **Pathway/Functional Enrichment**  
   - Tool: clusterProfiler  
   - Input: DEG list  
   - Output: Enriched GO terms / KEGG pathways  

### Step 7–8: Tumor Microenvironment Characterization
7. **Immune Signature Analysis**  
   - Tool: GSVA (ssGSEA method)  
   - Input: Normalized count matrix (vst/rlog) + immune-related gene sets  
   - Output: Per-sample immune signature enrichment scores  
8. **Immune Cell Deconvolution & Cell State Inference**  
   - Tools: CIBERSORTx, MCP-counter, xCell, Ecotyper  
   - Inputs:  
     - CIBERSORTx / xCell: TPM matrix  
     - MCP-counter: Raw count matrix  
     - Ecotyper: TPM or normalized count matrix (TPM recommended)  
   - Outputs:  
     - Immune cell proportions per sample  
     - Cell states and ecotypes (Ecotyper)  

### Final Integration
- Integrate DGEA and immune profiling results  
- Methods: correlation analysis, heatmaps, clustering  
- Final Output: comprehensive immunogram (TME characterization report)

> Each step provides input/output formats and usage notes in `docs/` or script headers.

## Licensing
- **Code**: MIT License (`LICENSE`)
- **Figures/Docs** (if needed): may use CC BY 4.0 (to be specified under `docs/`)

## Citation
> If you use this repository, please cite:  
> *Kobayashi Y., et al., “Bulk RNA-Seq-Based Deconvolution and Functional Profiling of the Tumor Immune Landscape: A Step-by-Step Analytical Guide” 2025.*  
(To be updated with publication details)

## Contact
- Maintainer: Yukari Kobayashi  
- Feedback via Issues/PR is welcome
