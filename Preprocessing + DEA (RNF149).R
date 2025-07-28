# Thesis pipeline step 2: Preprocessing and classifying samples
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
# Loading data:
data <- readRDS("TCGA_PAAD_RNAseq_STAR_Counts.rds")

expr_matrix <- assay(data)

head(rownames(expr_matrix))

# Removing versions:
noversion_gene_ids <- sub("\\..*", "", rownames(expr_matrix))

rownames(expr_matrix) <- noversion_gene_ids

# Checking for duplicates
duplicated_genes <- any(duplicated(noversion_gene_ids))

dup_gene_counts <- table(noversion_gene_ids)
dup_gene_counts <- dup_gene_counts[dup_gene_counts > 1]
dup_gene_counts

# Aggregating duplicates
expr_aggregated <- aggregate(expr_matrix, by = list(GeneID = rownames(expr_matrix)), FUN = mean)

rownames(expr_aggregated) <- expr_aggregated$GeneID
expr_aggregated <- expr_aggregated[, -1]

head(expr_aggregated)

# Save raw counts before log transformation [Crucial for Differential Expression Analysis (DESeq2 requirement)]
saveRDS(expr_aggregated, "TCGA_PAAD_RawCounts_Aggregated.rds")

# Log transformation
expr_transformed <- log2(expr_aggregated + 1)

head(expr_transformed)
str(expr_transformed)

expr_log2 <- data.matrix(expr_transformed)

# Find median RNF149 expression level
RNF149_expr <- expr_log2["ENSG00000163162", ]
RNF149_median <- median(RNF149_expr, na.rm = TRUE)

# Classifying samples according to median RNF149 expression
RNF149_group <- ifelse(RNF149_expr > RNF149_median, "High", "Low")

# Adding a new column to the metadata
colData(data)$RNF149_GROUP <- RNF149_group

head(colData(data))

table(colData(data)$RNF149_GROUP)

saveRDS(data, "TCGA_PAAD_Classification_Column_RNF149.rds")

# Thesis pipeline step 3: Differential Expression Analysis
# Load required packages
library(DESeq2)

# Load preprocessed data
data_2 <- readRDS("TCGA_PAAD_Classification_Column_RNF149.rds")

# Load the saved raw counts matrix (aggregated, integer counts)
raw_counts <- readRDS("TCGA_PAAD_RawCounts_Aggregated.rds")

# Round counts to integers (DESeq2 requirement)
raw_counts <- round(raw_counts)

# Extract sample metadata with RNF149_GROUP classification
coldata <- colData(data_2)

# Make sure RNF149_GROUP is a factor with correct levels
coldata$RNF149_GROUP <- factor(coldata$RNF149_GROUP, levels = c("Low", "High"))

# Create DESeqDataSet object with no-intercept model
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = coldata,
  design = ~ 0 + RNF149_GROUP  
)

# Run DESeq2 object
dds <- DESeq(dds)

# Extract results for High vs Low using the contrast argument (crucial for group-to-group comparisons)
# No reference group is set in this mode (both groups are explicit)
# The contrast function computes the difference between the means of the two groups 
res <- results(dds, contrast = c("RNF149_GROUP", "High", "Low"), alpha = 0.05)

# Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Ensembl ID to Gene Symbol conversion (for visualization later on)
install.packages("EnsDb.Hsapiens.v79_2.99.0.tar.gz", repos = NULL, type = "source")
library(EnsDb.Hsapiens.v79)
library(ensembldb)

# Assign the database object to a variable
edb <- EnsDb.Hsapiens.v79

# Convert DESeqResults to data frame for annotation only
# DESeq2 package update now prevents merging DESeqResults objects directly without data frames :(
res_ordered_df <- as.data.frame(res_ordered)

# Double check that the version numbers from Ensembl IDs are removed for annotation
res_ordered_df$ensembl <- sub("\\..*", "", rownames(res_ordered_df))

# Get mapping from EnsDb
mapping <- ensembldb::select(
  edb,
  keys = res_ordered_df$ensembl,
  keytype = "GENEID",
  columns = c("GENEID", "SYMBOL")
)

# Merge annotation with results
res_annotated <- merge(
  res_ordered_df, mapping,
  by.x = "ensembl", by.y = "GENEID",
  all.x = TRUE, sort = FALSE
)

# Calculate mapping loss for all genes
num_total <- nrow(res_annotated)
num_mapped <- sum(!is.na(res_annotated$SYMBOL) & res_annotated$SYMBOL != "")
num_unmapped <- num_total - num_mapped
percent_loss <- round(100 * num_unmapped / num_total, 2)
cat(sprintf("Total genes: %d\nMapped to gene symbol: %d\nUnmapped: %d (%.2f%% loss)\n",
            num_total, num_mapped, num_unmapped, percent_loss))

# Filter significant DEGs 
significant_degs_annot <- subset(res_annotated, padj < 0.05 & abs(log2FoldChange) > 1)

significant_degs_annot$Regulation <- ifelse(significant_degs_annot$log2FoldChange > 0, "Upregulated", "Downregulated")

# Calculate mapping loss for significant DEGs
num_sig_total <- nrow(significant_degs_annot)
num_sig_mapped <- sum(!is.na(significant_degs_annot$SYMBOL) & significant_degs_annot$SYMBOL != "")
num_sig_unmapped <- num_sig_total - num_sig_mapped
percent_sig_loss <- round(100 * num_sig_unmapped / num_sig_total, 2)
cat(sprintf("Significant DEGs: %d\nMapped to gene symbol: %d\nUnmapped: %d (%.2f%% loss)\n",
            num_sig_total, num_sig_mapped, num_sig_unmapped, percent_sig_loss))

# Save DEGs to a CSV file
write.csv(significant_degs_annot, "RNF149_Associated_DEGs_Annotated.csv", row.names = FALSE)

# Thesis pipeline step 3.1: Visualization
# Create the volcano plot (EnhancedVolcano)
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

# Create the Volcano plot using the DEA results object (annotated)
# Label using gene symbols
EnhancedVolcano(
  res_annotated,
  lab = res_annotated$SYMBOL,                  
  x = 'log2FoldChange',
  y = 'padj',
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'adjusted p-value'),
  title = 'Volcano Plot: RNF149 High vs Low',
  subtitle = 'TCGA PAAD RNA-seq',
  pCutoff = 0.05,                              # FDR threshold
  FCcutoff = 1,                                # log2FC threshold
  pointSize = 2.0,
  labSize = 3.5,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  legendLabels = c('NS','Log2FC','FDR','FDR & Log2FC'),
  legendPosition = 'right',
  drawConnectors = TRUE,                       
  widthConnectors = 0.5,
  max.overlaps = 30                            
)

# Save the Volcano plot to file 
# Use a graphic device (PDF/PNG) as ggsave does not work directly with EnhancedVolcano (not a ggplot object)
png("RNF149_VolcanoPlot.png", width=1200, height=900, res=150)

EnhancedVolcano(
  res_annotated,
  lab = res_annotated$SYMBOL,
  x = 'log2FoldChange',
  y = 'padj',
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'adjusted p-value'),
  title = 'Volcano Plot: RNF149 High vs Low',
  subtitle = 'TCGA PAAD RNA-seq',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 3.5,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  legendLabels = c('NS','Log2FC','FDR','FDR & Log2FC'),
  legendPosition = 'right',
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = 30
)
dev.off()

# Create the MA plot  using the DEA results object (annotated)
library(ggplot2)

# Colour genes based on the FDR threshold
res_annotated$Significant <- ifelse(res_annotated$padj < 0.05 & !is.na(res_annotated$padj), "FDR < 0.05", "NS")

# Label only the most significant DEGs to avoid overlapping
top_genes <- res_annotated[res_annotated$padj < 0.05 & abs(res_annotated$log2FoldChange) > 1 & !is.na(res_annotated$SYMBOL), ]
top_genes <- top_genes[order(top_genes$padj), ]
top_genes <- head(top_genes, 20)  

# Basic MA plot
ma_plot <- ggplot(res_annotated, aes(x = log10(baseMean + 1), y = log2FoldChange, color = Significant)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("FDR < 0.05" = "red", "NS" = "grey60")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "MA Plot: RNF149 High vs Low",
    x = "log10 Mean Expression",
    y = "Log2 Fold Change"
  ) +
  theme_bw() +
  theme(legend.position = "top")

# Add gene symbol labels for top genes
library(ggrepel)
ma_plot <- ma_plot + geom_text_repel(
  data = top_genes,
  aes(label = SYMBOL),
  size = 3,
  color = "black",
  box.padding = 0.3,
  max.overlaps = 20
)
print(ma_plot)

# Save the MA plot
ggsave("RNF149_MA_plot.png", plot = ma_plot, width = 8, height = 6, dpi = 150)

# Create the heatmap 
# Apply VST transformation for better visualization
vsd <- vst(dds, blind = FALSE)
vst_expr <- assay(vsd)

# Select top N significant DEGs by adjusted p-value
N <- 40

# Remove genes with NA or empty SYMBOL
top_degs_clean <- significant_degs_annot[
  !is.na(significant_degs_annot$SYMBOL) & significant_degs_annot$SYMBOL != "", ]

# Keep unique gene symbols (to avoid duplicates in heatmap rownames)
top_degs_clean <- top_degs_clean[!duplicated(top_degs_clean$SYMBOL), ]

# Order by significance
top_degs_clean <- top_degs_clean[order(top_degs_clean$padj), ]

# Select top N
top_degs <- head(top_degs_clean, N)
top_symbols <- top_degs$SYMBOL
top_ensembl <- top_degs$ensembl

# Filter out genes with very low variance across samples
gene_variances <- apply(vst_expr[top_ensembl, , drop = FALSE], 1, stats::var)
high_var_genes <- gene_variances > 0.1
top_ensembl <- top_ensembl[high_var_genes]
top_symbols <- top_symbols[high_var_genes]

# Label using the gene symbol
heatmap_mat <- vst_expr[top_ensembl, ]
rownames(heatmap_mat) <- top_symbols

# Create a dataframe for sample annotation (High, Low)
# Previously assigned to the metadata
annotation_col <- data.frame(
  RNF149_GROUP = colData(data)$RNF149_GROUP
)
rownames(annotation_col) <- colnames(heatmap_mat)

# Draw the heatmap
# scale = 'row' (z-score normalization) removes absolute expression differences
library(grid)
library(pheatmap)
heatmap_obj <- pheatmap(
  heatmap_mat,
  annotation_col = annotation_col,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 5,
  cellheight = 5,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top 40 DEGs Heatmap (RNF149 High vs Low)"
)

# Save the heatmap to file
save_pheatmap_png <- function(pheatmap_obj, filename, gene_count, width=1200, res=150) {
  height <- max(600, gene_count * 25)  
  png(filename, width=width, height=height, res=res)
  grid.newpage()
  grid.draw(pheatmap_obj$gtable)
  dev.off()
}

# Save to PNG
save_pheatmap_png(heatmap_obj, "RNF149_Top40_DEGs_Heatmap.png", gene_count = nrow(heatmap_mat))

