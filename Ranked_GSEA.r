# Extended Fold‐change GSEA: GO (BP/MF/CC), KEGG, WikiPathways, BioCarta

dependencies <- c(
  "BiocManager","dplyr","tibble","ggplot2",
  "clusterProfiler","org.Hs.eg.db","msigdbr","enrichplot"
)
for (pkg in dependencies) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    if (pkg %in% c("clusterProfiler","org.Hs.eg.db","msigdbr","enrichplot")) {
      BiocManager::install(pkg, ask=FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

library(dplyr)
library(tibble)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)

# 1) Read your DEG table
degs <- read.csv(
  "C:/Users/Nourine/Desktop/MBOAT7/MBOAT7_Associated_DEGs(No intercept)_Annotated(ensembldb).csv",
  stringsAsFactors = FALSE
) %>%
  mutate(
    SYMBOL = .[[8]],         # copy the 8th column into SYMBOL
    log2FC = log2FoldChange  
  ) %>%
  filter(!is.na(SYMBOL), !is.na(log2FC))

# 2) Map SYMBOL → ENTREZID, drop unmapped
degs$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys      = degs$SYMBOL,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)
degs <- dplyr::filter(degs, !is.na(ENTREZID))

# 3) Build a named, sorted fold‐change vector
rawList <- degs$log2FC
names(rawList) <- degs$ENTREZID
collapsed <- tapply(rawList, names(rawList), mean)
geneList  <- sort(as.numeric(collapsed), decreasing = TRUE)
names(geneList) <- names(collapsed)

# A) GSEA: GO and KEGG via clusterProfiler
gsea_bp   <- gseGO(
  geneList, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
  ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE
)
gsea_mf   <- gseGO(
  geneList, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
  ont="MF", pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE
)
gsea_cc   <- gseGO(
  geneList, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
  ont="CC", pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE
)
gsea_kegg <- gseKEGG(
  geneList, organism="hsa",
  pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE
)

# B) GSEA: WikiPathways & BioCarta via msigdbr + clusterProfiler::GSEA
# Use updated msigdbr args: collection, subcollection
wp_sets <- msigdbr(
  species       = "Homo sapiens",
  collection    = "C2",
  subcollection = "CP:WIKIPATHWAYS"
) %>% dplyr::select(gs_name, ncbi_gene)
gsea_wp <- GSEA(
  geneList      = geneList,
  TERM2GENE     = wp_sets,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

bc_sets <- msigdbr(
  species       = "Homo sapiens",
  collection    = "C2",
  subcollection = "CP:BIOCARTA"
) %>% dplyr::select(gs_name, ncbi_gene)
gsea_bc <- GSEA(
  geneList      = geneList,
  TERM2GENE     = bc_sets,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

# C) Plotting helper: manual ggplot2 barplot of top‐10 GSEA results
plot_top10 <- function(gsea_res, title) {
  df <- as.data.frame(gsea_res@result) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10) %>%
    mutate(Description = factor(Description, levels = rev(Description)))
  ggplot(df, aes(x = NES, y = Description, fill = p.adjust)) +
    geom_col() +
    scale_fill_viridis_c(option = "plasma", direction = -1) +
    labs(title = title, x = "Normalized Enrichment Score", y = NULL, fill = "adj. p-value") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 10))
}

# D) Generate & print four barplots
print(plot_top10(gsea_bp,   "GO BP"))
print(plot_top10(gsea_mf,   "GO MF"))
print(plot_top10(gsea_cc,   "GO CC"))
print(plot_top10(gsea_kegg, "KEGG"))
print(plot_top10(gsea_wp,   "WikiPathways"))
print(plot_top10(gsea_bc,   "BioCarta"))

# E) Save full GSEA result tables
#write.csv(as.data.frame(gsea_bp@result),   "GSEA_GO_BP.csv",         row.names = FALSE)
#write.csv(as.data.frame(gsea_mf@result),   "GSEA_GO_MF.csv",         row.names = FALSE)
#write.csv(as.data.frame(gsea_cc@result),   "GSEA_GO_CC.csv",         row.names = FALSE)
#write.csv(as.data.frame(gsea_kegg@result), "GSEA_KEGG.csv",          row.names = FALSE)
#write.csv(as.data.frame(gsea_wp@result),   "GSEA_WikiPathways.csv",  row.names = FALSE)
#write.csv(as.data.frame(gsea_bc@result),   "GSEA_BioCarta.csv",      row.names = FALSE)

#########################
# Test code

# Helper to convert a result@result data.frame's core_enrichment from ENTREZ to SYMBOL
convert_entrez_to_symbol <- function(gsea_df) {
  gsea_df$core_enrichment <- sapply(gsea_df$core_enrichment, function(ent_str) {
    entrez_ids <- strsplit(ent_str, "/")[[1]]
    syms <- mapIds(
      org.Hs.eg.db,
      keys      = entrez_ids,
      column    = "SYMBOL",
      keytype   = "ENTREZID",
      multiVals = "first"
    )
    # Remove any NAs and collapse back into one string
    paste(na.omit(syms), collapse = "/")
  })
  gsea_df
}

# D) Generate & print four barplots (unchanged)
print(plot_top10(gsea_bp,   "GO BP (Top 10)"))
print(plot_top10(gsea_mf,   "GO MF (Top 10)"))
print(plot_top10(gsea_cc,   "GO CC (Top 10)"))
print(plot_top10(gsea_kegg, "KEGG (Top 10)"))
print(plot_top10(gsea_wp,   "WikiPathways (Top 10)"))
print(plot_top10(gsea_bc,   "BioCarta (Top 10)"))

# E) Save full GSEA result tables with SYMBOLs in core_enrichment
# First extract and convert each result
res_bp   <- convert_entrez_to_symbol(as.data.frame(gsea_bp@result))
res_mf   <- convert_entrez_to_symbol(as.data.frame(gsea_mf@result))
res_cc   <- convert_entrez_to_symbol(as.data.frame(gsea_cc@result))
res_kegg <- convert_entrez_to_symbol(as.data.frame(gsea_kegg@result))
res_wp   <- convert_entrez_to_symbol(as.data.frame(gsea_wp@result))
res_bc   <- convert_entrez_to_symbol(as.data.frame(gsea_bc@result))

# Then write them out
write.csv(res_bp,   "mboat7GSEA_GO_BP_symbols.csv",        row.names = FALSE)
write.csv(res_mf,   "mboat7GSEA_GO_MF_symbols.csv",        row.names = FALSE)
write.csv(res_cc,   "mboat7GSEA_GO_CC_symbols.csv",        row.names = FALSE)
write.csv(res_kegg, "mboat7GSEA_KEGG_symbols.csv",         row.names = FALSE)
write.csv(res_wp,   "mboat7GSEA_WikiPathways_symbols.csv", row.names = FALSE)
write.csv(res_bc,   "mboat7GSEA_BioCarta_symbols.csv",     row.names = FALSE)
