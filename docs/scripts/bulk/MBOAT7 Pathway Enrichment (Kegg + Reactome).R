###############################################################################
# GSEA from Ensembl-based DEGs using log2FoldChange (KEGG + Reactome)
# Improved: jitter to fix ties, eps=0 for p-values, parallelization kept
###############################################################################

## 1 ── Load libraries ----
needed_pkgs <- c("BiocManager",
                 "clusterProfiler", "ReactomePA", "org.Hs.eg.db",
                 "enrichplot", "data.table", "AnnotationDbi", "biomaRt")

for (pkg in needed_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "BiocManager") install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(enrichplot)
library(data.table)
library(biomaRt)
library(clusterProfiler)
library(enrichplot)


## 2 ── Load DEG file ----
file_path <- "D:/Downloads/MBOAT7_Associated_DEGs_Annotated.csv"
deg <- fread(file_path)

## 3 ── Map Ensembl IDs to Entrez IDs ----
ensembl_ids <- unique(deg$ensembl)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = ensembl_ids,
                 mart = mart)

setnames(mapping, old = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
         new = c("ensembl", "ENTREZID", "SYMBOL"))

deg_mapped <- merge(deg, mapping, by = "ensembl")
deg_mapped <- deg_mapped[!is.na(ENTREZID)]

## 4 ── Prepare log2FoldChange safely ----
lfc_col <- grep("log2FoldChange", names(deg_mapped), value = TRUE, ignore.case = TRUE)

if (length(lfc_col) == 0) {
  stop("❌ 'log2FoldChange' column not found.")
} else if (length(lfc_col) > 1) {
  cat("⚠ Multiple 'log2FoldChange' columns found:\n")
  print(lfc_col)
  lfc_col <- lfc_col[1]
}
setnames(deg_mapped, lfc_col, "log2FoldChange")

# Add jitter to fix GSEA warning about tied stats
deg_mapped$log2FoldChange <- deg_mapped$log2FoldChange + rnorm(nrow(deg_mapped), 0, 1e-6)

# Deduplicate by most significant fold change per Entrez
deg_mapped <- deg_mapped[order(-abs(log2FoldChange))]
deg_mapped <- deg_mapped[!duplicated(deg_mapped$ENTREZID)]

## 5 ── Create ranked gene list ----
geneList <- deg_mapped$log2FoldChange
names(geneList) <- deg_mapped$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

## 6 ── Run GSEA for KEGG ----
# Note: Parallelization is kept enabled; serialize warnings can be ignored
gsea_kegg <- gseKEGG(geneList     = geneList,
                     organism     = "hsa",
                     keyType      = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     eps          = 0)

## 7 ── Run GSEA for Reactome ----
gsea_react <- gsePathway(geneList     = geneList,
                         organism     = "human",
                         pvalueCutoff = 0.05,
                         minGSSize    = 10,
                         maxGSSize    = 500,
                         eps          = 0)


## 8 ── Save RESULTS

fwrite(as.data.table(gsea_kegg),  "GSEA_mboat7_KEGG_results.csv")
fwrite(as.data.table(gsea_react), "GSEA_mboat7_Reactome_results.csv")

## 9 ── Show core enrichment genes for KEGG ----
core_entrez <- unique(unlist(strsplit(gsea_kegg@result$core_enrichment, "/")))
deg_core_overlap <- deg_mapped[ENTREZID %in% core_entrez]
View(deg_core_overlap)

