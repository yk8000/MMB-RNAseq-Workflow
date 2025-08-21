## Tools
EnhancedVolcano (R package)

## Inputs
DEG_results_raw.csv : DESeq2 results from Step5-1 (rows = genes, columns include log2FoldChange, padj)

## Outputs
volcano_pre_vs_post.png : Volcano plot visualizing DEGs (post vs pre)

## Key Options
pCutoff = 0.05 : adjusted p-value threshold
FCcutoff = 1 : absolute log2 fold change threshold
selectLab : top 10 significant up- and down-regulated genes labeled
col : customized color scheme (grey50, forestgreen, royalblue, red2)
titleLabSize / axisLabSize / labSize : adjust font sizes for clarity
drawConnectors = TRUE : link gene labels to points for readability

## Notes
This plot directly visualizes results from Step5-1 (DESeq2 DGE analysis).
By default, labels are Ensembl IDs. If org.Hs.eg.db is available, they are mapped to gene symbols.
Significant DEGs are defined as padj < 0.05 and |log2FC| > 1.
Among significant DEGs, the top 10 genes (per direction) are labeled based on smallest adjusted p-values.
Adjust thresholds (pCutoff, FCcutoff) depending on desired stringency.
