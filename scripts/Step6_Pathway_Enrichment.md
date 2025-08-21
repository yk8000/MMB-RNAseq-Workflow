# Step6. Pathway and Functional Enrichment Analysis with clusterProfiler
This script performs Over-Representation Analysis (ORA) using **clusterProfiler**.  
It takes DEG results from Step5-1 (`DEG_results_raw.csv`), separates UP and DOWN genes,  
runs GO (BP/CC/MF) and KEGG enrichment, and outputs CSV tables and dot plots.

```
# Rscript

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)

set.seed(1234)  # project convention

# ---------------------------
# 0) Parameters
# ---------------------------
# - alpha: significance threshold (adjusted p-value)
# - lfc:   absolute log2 fold change cutoff
# - infile: DEG results file from Step5-1
# - species_kegg: KEGG organism code
alpha <- 0.05
lfc   <- 1
species_kegg <- "hsa"
infile <- "DEG_results_raw.csv"

topn_go_bp  <- 10
topn_go_etc <- 10
topn_kegg   <- 10

# ---------------------------
# 1) Load DE results
# ---------------------------
# - Input must contain 'padj' and 'log2FoldChange'
# - Replace missing padj values with 1 (to treat as non-significant)
res <- read.csv(infile, row.names = 1, check.names = FALSE) |> as.data.frame()
if (!all(c("padj","log2FoldChange") %in% colnames(res))) {
  stop("Input must contain columns 'padj' and 'log2FoldChange'.")
}
res$padj[is.na(res$padj)] <- 1

# ---------------------------
# 2) Define gene sets (UP/DOWN)
# ---------------------------
# - UP: padj < alpha & log2FC > lfc
# - DOWN: padj < alpha & log2FC <= -lfc (inclusive boundary)
up_ids_raw   <- rownames(res)[res$padj < alpha & res$log2FoldChange >  lfc]
down_ids_raw <- rownames(res)[res$padj < alpha & res$log2FoldChange <= -lfc]

# ---------------------------
# 3) Build GO universes and sets (ENSEMBL)
# ---------------------------
# - strip version suffix from Ensembl IDs (e.g., ENSG00000141510.12 → ENSG00000141510)
# - create GO universe and UP/DOWN sets
strip_version <- function(x) sub("\\..*$", "", x)

ens_universe <- unique(strip_version(rownames(res)))
ens_up       <- unique(strip_version(up_ids_raw))
ens_down     <- unique(strip_version(down_ids_raw))

# ---------------------------
# 4) Map Ensembl -> ENTREZ for KEGG
# ---------------------------
# - Use clusterProfiler::bitr to convert Ensembl to ENTREZ + SYMBOL
# - Drop unmapped IDs
safe_df <- function(x) if (is.null(x)) data.frame() else x

map_entrez <- tryCatch(
  bitr(unique(ens_universe), fromType = "ENSEMBL",
       toType = c("ENTREZID","SYMBOL"), OrgDb = org.Hs.eg.db),
  error = function(e) data.frame()
) |> safe_df()

to_entrez <- function(ens_vec, map_df) {
  if (!length(ens_vec)) return(tibble(ENTREZID = character(), SYMBOL = character()))
  tibble(ENSEMBL = ens_vec) |>
    left_join(map_df, by = c("ENSEMBL" = "ENSEMBL")) |>
    filter(!is.na(ENTREZID)) |>
    distinct(ENTREZID, .keep_all = TRUE)
}

universe_tbl <- to_entrez(ens_universe, map_entrez)
up_tbl       <- to_entrez(ens_up,       map_entrez)
down_tbl     <- to_entrez(ens_down,     map_entrez)

gene_universe_kegg <- unique(universe_tbl$ENTREZID)
genes_up_kegg      <- unique(up_tbl$ENTREZID)
genes_down_kegg    <- unique(down_tbl$ENTREZID)

message("GO universe (ENSEMBL): ", length(ens_universe))
message("KEGG universe (ENTREZ): ", length(gene_universe_kegg))
message("UP: ENSEMBL=", length(ens_up), " | ENTREZ=", length(genes_up_kegg))
message("DOWN: ENSEMBL=", length(ens_down), " | ENTREZ=", length(genes_down_kegg))

# ---------------------------
# 5) Output utilities
# ---------------------------
# - mkout: create output dirs (UP, DOWN)
# - save_tbl: save tables as CSV
# - save_empty_set: output empty CSVs when no DEGs
# - plot_dot: draw dot plots of enrichment results
# - simplify_bp: reduce redundancy in GO BP terms
mkout <- function(dir) if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
mkout("UP"); mkout("DOWN")

save_tbl <- function(x, path) {
  if (is.null(x)) x <- data.frame()
  write.csv(x, path, row.names = FALSE)
}
save_empty_set <- function(side) {
  save_tbl(data.frame(), file.path(side, paste0("GO_BP_", side, ".csv")))
  save_tbl(data.frame(), file.path(side, paste0("GO_CC_", side, ".csv")))
  save_tbl(data.frame(), file.path(side, paste0("GO_MF_", side, ".csv")))
  save_tbl(data.frame(), file.path(side, paste0("KEGG_",  side, ".csv")))
}
plot_dot <- function(edo, outfile, topn = 20, title = NULL) {
  df <- tryCatch(as.data.frame(edo), error = function(e) data.frame())
  if (nrow(df) == 0) return(invisible())
  png(outfile, width = 2000, height = 1800, res = 300)
  print(dotplot(edo, showCategory = topn, title = title))
  dev.off()
}
simplify_bp <- function(edo, cutoff = 0.6) {
  df <- tryCatch(as.data.frame(edo), error = function(e) data.frame())
  if (nrow(df) == 0) return(edo)
  tryCatch(simplify(edo, cutoff = cutoff, by = "p.adjust", select_fun = min),
           error = function(e) edo)
}

# ---------------------------
# 6) Enrichment wrappers
# ---------------------------
# - run_enrich_go_ens: GO ORA with Ensembl IDs (results readable in SYMBOLs)
# - run_enrich_kegg: KEGG ORA with ENTREZ IDs
run_enrich_go_ens <- function(ens_gene_ids, ens_univ_ids, ont = "BP") {
  enrichGO(gene          = ens_gene_ids,
           universe      = ens_univ_ids,
           OrgDb         = org.Hs.eg.db,
           keyType       = "ENSEMBL",
           ont           = ont,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.2,
           readable      = TRUE)
}
run_enrich_kegg <- function(entrez_gene_ids, entrez_univ_ids, organism = "hsa") {
  enrichKEGG(gene          = entrez_gene_ids,
             universe      = entrez_univ_ids,
             organism      = organism,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.2)
}

# ---------------------------
# 7) Run for one side
# ---------------------------
# - Save DEG lists (Ensembl + mapped ENTREZ)
# - Run GO and KEGG enrichment
# - Save CSVs and dot plots
# - If no DEGs, output empty files
one_side <- function(side, ens_ids_side, entrez_tbl_side) {
  outdir <- side
  # Save gene lists
  write.csv(data.frame(ENSEMBL = ens_ids_side),
            file.path(outdir, paste0(tolower(side), "_genes_used_GO_ENSEMBL.csv")),
            row.names = FALSE)
  save_tbl(entrez_tbl_side, file.path(outdir, paste0(tolower(side), "_genes_mapped_ENTREZ.csv")))

  if (!length(ens_ids_side) && nrow(entrez_tbl_side) == 0) {
    message(side, ": no genes; writing empty CSVs.")
    save_empty_set(side)
    return(invisible())
  }

  # --- GO (ENSEMBL) ---
  if (length(ens_ids_side)) {
    ego_bp <- run_enrich_go_ens(ens_ids_side, ens_universe, "BP") |> simplify_bp(0.6)
    ego_cc <- run_enrich_go_ens(ens_ids_side, ens_universe, "CC")
    ego_mf <- run_enrich_go_ens(ens_ids_side, ens_universe, "MF")

    save_tbl(as.data.frame(ego_bp), file.path(outdir, paste0("GO_BP_", side, ".csv")))
    save_tbl(as.data.frame(ego_cc), file.path(outdir, paste0("GO_CC_", side, ".csv")))
    save_tbl(as.data.frame(ego_mf), file.path(outdir, paste0("GO_MF_", side, ".csv")))

    plot_dot(ego_bp, file.path(outdir, paste0("GO_BP_", side, "_dotplot.png")),
             topn = topn_go_bp, title = paste0("GO BP (", side, ")"))
    plot_dot(ego_cc, file.path(outdir, paste0("GO_CC_", side, "_dotplot.png")),
             topn = topn_go_etc, title = paste0("GO CC (", side, ")"))
    plot_dot(ego_mf, file.path(outdir, paste0("GO_MF_", side, "_dotplot.png")),
             topn = topn_go_etc, title = paste0("GO MF (", side, ")"))
  } else {
    save_tbl(data.frame(), file.path(outdir, paste0("GO_BP_", side, ".csv")))
    save_tbl(data.frame(), file.path(outdir, paste0("GO_CC_", side, ".csv")))
    save_tbl(data.frame(), file.path(outdir, paste0("GO_MF_", side, ".csv")))
  }

  # --- KEGG (ENTREZ) ---
  if (nrow(entrez_tbl_side)) {
    genes_entrez <- unique(entrez_tbl_side$ENTREZID)
    ekegg <- run_enrich_kegg(genes_entrez, gene_universe_kegg, species_kegg)
    save_tbl(as.data.frame(ekegg), file.path(outdir, paste0("KEGG_", side, ".csv")))
    plot_dot(ekegg, file.path(outdir, paste0("KEGG_", side, "_dotplot.png")),
             topn = topn_kegg, title = paste0("KEGG (", side, ")"))
  } else {
    save_tbl(data.frame(), file.path(outdir, paste0("KEGG_", side, ".csv")))
  }
}

# ---------------------------
# 8) Execute (UP and DOWN)
# ---------------------------
# - Run enrichment for UP and DOWN DEGs separately
one_side("UP",   ens_up,   up_tbl)
one_side("DOWN", ens_down, down_tbl)
