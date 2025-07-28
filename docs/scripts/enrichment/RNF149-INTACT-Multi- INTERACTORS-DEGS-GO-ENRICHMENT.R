# 📦 Step 0: Install & load required packages
required_packages <- c("httr", "jsonlite", "dplyr", "clusterProfiler",
                       "org.Hs.eg.db", "ReactomePA", "OmnipathR", "stringr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if(pkg == "OmnipathR") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("OmnipathR")
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# ------------------------
# 🧬 Step 1: Query IntAct (PSICQUIC) for RNF149 interactors
cat("🔎 Querying IntAct for RNF149 interactors...\n")
query_gene <- "RNF149"
psicquic_url <- paste0(
  "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/",
  query_gene,
  "?format=tab25"
)

response <- tryCatch(httr::GET(psicquic_url), error = function(e) NULL)
if (is.null(response) || httr::status_code(response) != 200) {
  stop("❌ Failed to retrieve data from IntAct/PSICQUIC.")
}

interactions_text <- httr::content(response, as = "text", encoding = "UTF-8")
interactions_df <- read.delim(text = interactions_text, header = FALSE, sep = "\t", quote = "")

extract_uniprot_id <- function(x) {
  stringr::str_extract(x, "uniprotkb:[^|]+") %>%
    gsub("uniprotkb:", "", .)
}

uniprot_ids_1 <- sapply(interactions_df$V1, extract_uniprot_id)
uniprot_ids_2 <- sapply(interactions_df$V2, extract_uniprot_id)
intact_uniprot <- unique(c(uniprot_ids_1, uniprot_ids_2))
intact_uniprot <- intact_uniprot[!is.na(intact_uniprot)]

cat("✅ IntAct interactors found:", length(intact_uniprot), "\n")

# Map UniProt IDs to gene SYMBOLs
intact_symbols <- bitr(intact_uniprot, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
intact_symbols <- unique(intact_symbols$SYMBOL)

# ------------------------
# 🧬 Step 2: Query OmniPath (includes MINT, HPRD, DIP, and others) for RNF149 interactors
cat("🔎 Querying OmniPath for RNF149 interactors (MINT, HPRD, DIP included)...\n")
omnipath_interactions <- import_omnipath_interactions()

# Filter interactions where RNF149 is source or target
omnipath_subset <- omnipath_interactions[
  omnipath_interactions$source_genesymbol == query_gene | omnipath_interactions$target_genesymbol == query_gene, ]

omnipath_interactors <- unique(c(omnipath_subset$source_genesymbol, omnipath_subset$target_genesymbol))
omnipath_interactors <- setdiff(omnipath_interactors, query_gene)
cat("✅ OmniPath interactors found:", length(omnipath_interactors), "\n")

# ------------------------
# 🧮 Step 3: Combine all interactors (union)
all_interactors <- unique(c(intact_symbols, omnipath_interactors))
cat("🔗 Total unique interactors combined:", length(all_interactors), "\n")

# ------------------------
# 🧾 Step 4: Load your DEGs list (replace with your actual data frame)
# Example placeholder:
# significant_degs_with_rnf149 <- data.frame(SYMBOL = c("TTYH1", "TMEM74", "RET", "SFTPC", "DLK1", "RNF149"))
degs_symbols <- significant_degs_with_rnf149$SYMBOL

# Find overlap between DEGs and interactors
overlapping_genes <- intersect(degs_symbols, all_interactors)
cat("✅ Overlapping genes with DEGs:", length(overlapping_genes), "\n")
print(overlapping_genes)

# ------------------------
# 🧬 Step 5: Enrichment analysis on overlapping genes
if(length(overlapping_genes) >= 3) {
  entrez_map <- bitr(overlapping_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- entrez_map$ENTREZID
  
  # GO Biological Process enrichment
  go_results <- enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  # Reactome pathway enrichment
  reactome_results <- enrichPathway(
    gene         = entrez_ids,
    organism     = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable     = TRUE
  )
  
  # ------------------------
  # 📊 Display top enrichment results
  cat("\n📊 Top GO Biological Process enrichment terms:\n")
  print(head(go_results@result, 10))
  
  cat("\n📊 Top Reactome pathway enrichment terms:\n")
  print(head(reactome_results@result, 10))
  
  # ------------------------
  # 📉 Plot enrichment barplots
  if(nrow(go_results@result) > 0) {
    barplot(go_results, showCategory = 10, title = "GO Biological Process Enrichment")
  }
  if(nrow(reactome_results@result) > 0) {
    barplot(reactome_results, showCategory = 10, title = "Reactome Pathway Enrichment")
  }
  
  # ------------------------
  # 💾 Save results
  output_dir <- "D:/Downloads/rnf149_MultiDB_output"
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  saveRDS(list(
    interactors       = all_interactors,
    overlapping_genes = overlapping_genes,
    go_results        = go_results,
    reactome_results  = reactome_results
  ), file = file.path(output_dir, "rnf149_MultiDB_results.rds"))
  
  write.csv(data.frame(SYMBOL = all_interactors),
            file = file.path(output_dir, "all_interactors.csv"), row.names = FALSE)
  write.csv(data.frame(SYMBOL = overlapping_genes),
            file = file.path(output_dir, "DEG_interactors.csv"), row.names = FALSE)
  write.csv(go_results@result,
            file = file.path(output_dir, "GO_enrichment.csv"), row.names = FALSE)
  write.csv(reactome_results@result,
            file = file.path(output_dir, "Reactome_enrichment.csv"), row.names = FALSE)
  
  cat("✅ All results saved to:", output_dir, "\n")
} else {
  cat("⚠️ Not enough overlapping genes (less than 3) for enrichment analysis.\n")
}
