# Load Required Library
library(Seurat)
library(DoubletFinder)
library(parallel) # detectCores()
library(ggplot2)
library(dplyr)
library(harmony)
library(SingleR)
library(SingleCellExperiment)


# 0. Load scRNA-seq data ----------

# Step A: Define paths and extract .tar.gz files
raw_path <- "/Users/Sajda/Desktop/Grad/GSE155698_RAW"
tar_files <- list.files(raw_path, pattern = "\\.tar\\.gz$", full.names = TRUE)
extract_dir <- file.path(raw_path, "untarred")
dir.create(extract_dir, showWarnings = FALSE)

# Extract tar files
for (f in tar_files) {
  sample_name <- tools::file_path_sans_ext(basename(f))
  untar(f, exdir = file.path(extract_dir, sample_name))
}

# Step B: Identify tumor and normal directories
all_sample_dirs <- list.dirs(extract_dir, recursive = FALSE, full.names = TRUE)
normal_dirs <- grep("AdjNorm", all_sample_dirs, value = TRUE)
tumor_dirs  <- grep("PDAC", all_sample_dirs, value = TRUE)

# Step C: Function to load Seurat objects with proper sample labeling
load_group <- function(dirs, label) {
  objects <- list()
  for (d in dirs) {
    tar_name <- basename(d)
    sample_name <- gsub("\\.tar$", "", tar_name)
    sub_dirs <- list.dirs(d, recursive = TRUE)
    matrix_path <- sub_dirs[grepl("filtered_(feature|gene)_bc_matrix/?$", sub_dirs)] 
    h5_file <- list.files(d, pattern = "filtered_.*\\.h5$", recursive = TRUE, full.names = TRUE)
    
    # Read matrix
    if (length(matrix_path) == 1 && dir.exists(matrix_path)) {
      data <- Read10X(data.dir = matrix_path)
    } else if (length(h5_file) == 1) {
      data <- Read10X_h5(h5_file)
    } else {
      warning(paste("No valid matrix for:", d))
      next
    }
    
    # Create Seurat object
    obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200, project = sample_name)
    obj$group <- label
    obj$orig.ident <- sample_name  # Correct orig.ident
    objects[[sample_name]] <- obj
  }
  return(objects)
}

# Load all Seurat objects
tumor_objs <- load_group(tumor_dirs, "tumor")
normal_objs <- load_group(normal_dirs, "normal")

# Step D: Merge tumor and normal Seurat objects
if (length(tumor_objs) > 1) {
  tumor_merged <- merge(tumor_objs[[1]], y = tumor_objs[-1], add.cell.ids = names(tumor_objs))
} else if (length(tumor_objs) == 1) {
  tumor_merged <- tumor_objs[[1]]
} else {
  stop("No tumor samples loaded")
}

if (length(normal_objs) > 1) {
  normal_merged <- merge(normal_objs[[1]], y = normal_objs[-1], add.cell.ids = names(normal_objs))
} else if (length(normal_objs) == 1) {
  normal_merged <- normal_objs[[1]]
} else {
  stop("No normal samples loaded")
}

# Final merge
merged_seurat <- merge(normal_merged, y = tumor_merged)
# Optional check
print(table(merged_seurat$group, merged_seurat$orig.ident))

# 1. Quality Control  ---------

# MT%
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio

# Calculate the average for nCount_RNA
mean_nCount_RNA <- mean(merged_seurat@meta.data$nCount_RNA)

# Calculate the average for nFeature_RNA
mean_nFeature_RNA <- mean(merged_seurat@meta.data$nFeature_RNA)

# Calculate the average for mitoRatio
mean_mitoRatio <- mean(merged_seurat$mitoRatio)

# Display the averages
mean_nCount_RNA
mean_nFeature_RNA
mean_mitoRatio

summary(merged_seurat$nCount_RNA)

# Visualize QC metrics as a violin plot before and after filtering
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), pt.size=0,group.by = "group")

# Visualize feature-feature relation using scatterplot before and after filtering
options(repr.plot.width=12, repr.plot.height=8)

# Create the plots
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "mitoRatio", group.by = "group")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "group")

# Combine the plots
plot1 + plot2


# 2. Filtering  ---------
table(merged_seurat$orig.ident)
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & mitoRatio < 10)
table(merged_seurat$orig.ident)


# 3. Normalize data  ---------
merged_seurat <- NormalizeData(object = merged_seurat, normalization.method = "LogNormalize", 
                               scale.factor = 10000)
#str(merged_seurat)

# 4. Identify highly variable features  ---------
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_seurat), 10)

# 5. Scaling  ---------
all.genes <- rownames(merged_seurat)
merged_seurat <- ScaleData(merged_seurat, features = all.genes, vars.to.regress = "mitoRatio")
#merged_seurat <- ScaleData(merged_seurat, features = all.genes)

#str(merged_seurat)

# 6. Perform Linear dimensionality reduction  ---------
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))

# Save the Seurat object after doublet detection
#saveRDS(merged_seurat, file = "before_doublet.rds")
# Load the Seurat object after doublet detection
#merged_seurat <- readRDS("before_doublet.rds")

# 7. Doublet Detection  ---------

# set seed
set.seed(8) # for reproducibility

# work in parallel (Use available cores)
options(mc.cores = detectCores() - 1)

merged.split <- SplitObject(merged_seurat, split.by = "orig.ident")
# loop through samples to find doublets
for (i in 1:length(merged.split)) {
  # print the sample we are on
  print(paste0("Sample ",i)) 
  merged.split[[i]] <- JoinLayers(merged.split[[i]])
  # Pre-process seurat object with standard seurat workflow
  padc.sample <- NormalizeData(merged.split[[i]])
  padc.sample <- FindVariableFeatures(padc.sample)
  padc.sample <- ScaleData(padc.sample)
  padc.sample <- RunPCA(padc.sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdev <- padc.sample@reductions$pca@stdev
  var <- stdev^2
  
  EndVar = 0
  
  for(j in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:j])
    expvar <- numerator/total
    if(EndVar == 0){
      if(expvar > 0.9){
        EndVar <- EndVar + 1
        PCNum <- j
      }
    }
  }
  #Confirm #PC's determined explain > 90% of variance
  sum(var[1:PCNum])/ sum(var)
  
  # finish pre-processing
  padc.sample <- RunUMAP(padc.sample, dims = 1:PCNum)
  padc.sample <- FindNeighbors(object = padc.sample, dims = 1:PCNum)              
  padc.sample <- FindClusters(object = padc.sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(padc.sample, PCs = 1:PCNum, num.cores = detectCores() - 1)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- padc.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  n.cells <- ncol(padc.sample)
  doublet.rate <- 0.008 * (n.cells / 1000)
  nExp.poi <- round(n.cells * doublet.rate)
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  padc.sample <- doubletFinder(seu = padc.sample, 
                               PCs = 1:PCNum, 
                               pK = optimal.pk,
                               nExp = nExp.poi.adj)
  
  df_col <- grep("DF.classifications", colnames(padc.sample@meta.data), value = TRUE)
  padc.sample$doublet_finder <- padc.sample@meta.data[[df_col]]
  
  # subset and save
  merged.singlets <- subset(padc.sample, doublet_finder == "Singlet")
  merged.split[[i]] <- merged.singlets
  remove(merged.singlets)
}

# converge merged.split
merged.singlets <- merge(x = merged.split[[1]], 
                         y = merged.split[2:length(merged.split)],
                         project = "PADC scRNAseq")

# compare before and after doublet finder
cat("Before DoubletFinder:", ncol(merged_seurat), "cells\n")
cat("After DoubletFinder:", ncol(merged.singlets), "cells\n")
table(merged_seurat$group)
table(merged.singlets$group)

# 7.1. re-run the global workflow on merged_singlets ---------
merged.singlets <- NormalizeData(merged.singlets)
merged.singlets <- FindVariableFeatures(merged.singlets, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged.singlets)
merged.singlets <- ScaleData(merged.singlets, features = all.genes)
merged.singlets <- RunPCA(merged.singlets, features = VariableFeatures(object = merged.singlets))

#Determine Dimensions for 90% Variance
stdev <- merged.singlets@reductions$pca@stdev
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}
#Confirm #PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)

# determine dimensionality of the data
ElbowPlot(merged.singlets, ndims = 50)

# Save the Seurat object after doublet detection
#saveRDS(merged.singlets, file = "merged.singlets.rds")
# Load the Seurat object after doublet detection
#merged.singlets <- readRDS("merged.singlets.rds")

# 8. Batch Correction using Harmony  ---------

#Find Neighbors + Find CLusters (without harmony batch correction)
merged.singlets <- FindNeighbors(object = merged.singlets, dims = 1:PCNum)
merged.singlets <- FindClusters(object = merged.singlets, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot (without harmony batch correction)
merged.singlets <- RunUMAP(object = merged.singlets, dims = 1:PCNum)
DimPlot(object = merged.singlets, reduction = "umap", label = TRUE, pt.size = 0.5)
merged.singlets[["UMAP_Clusters"]] <- Idents(object = merged.singlets)

#Optionally annotate by group (e.g. tumor vs. normal)
DimPlot(merged.singlets, reduction = "umap",split.by = "group",label = TRUE, repel = TRUE, pt.size = 0.5)

#Find Neighbors + Find CLusters (with harmony batch correction)
# Run Harmony to remove batch effect from PCA space
merged.singlets <- RunHarmony(merged.singlets, group.by.vars = "orig.ident", plot_convergence = TRUE)

# Clustering on Harmony embeddings
merged_seurat.harmony <- FindNeighbors(merged.singlets, dims = 1:PCNum, reduction = "harmony")
merged_seurat.harmony <- FindClusters(merged.singlets, resolution = 1.2, reduction ="harmony")

# Run UMAP on Harmony embeddings
merged_seurat.harmony <- RunUMAP(merged_seurat.harmony, dims = 1:PCNum, reduction = "harmony")

# Plot UMAP
DimPlot(merged_seurat.harmony, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)
DimPlot(merged_seurat.harmony, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "orig.ident", ncol = 5)

# Store cluster identities
merged_seurat.harmony[["UMAP_Clusters"]] <- Idents(merged_seurat.harmony)

# Reset back to cluster IDs if needed
Idents(merged_seurat.harmony) <- "UMAP_Clusters"
levels(merged_seurat.harmony)

#Optionally annotate by group (e.g. tumor vs. normal)
DimPlot(merged_seurat.harmony, reduction = "umap",split.by = "group",label = TRUE, repel = TRUE, pt.size = 0.5)

# Save the Seurat object after doublet detection
#saveRDS(merged_seurat.harmony, file = "merged_seurat.harmony.rds")
# Load the Seurat object after doublet detection
#merged_seurat.harmony <- readRDS("merged_seurat.harmony.rds")

# 9. Cluster biomarkers  ---------
# Find all markers

merged_seurat.harmony <- JoinLayers(merged_seurat.harmony)
merged_seurat.harmony.markers <- FindAllMarkers(
  merged_seurat.harmony,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 1.5)


# 10. SingleR Annotation  ---------

# Load full PDAC reference
pdac_reference_seurat <- readRDS("~/Grad/pk_all.rds")

# Build SingleCellExperiment for reference
reference_expression <- GetAssayData(pdac_reference_seurat, assay = "RNA", layer  = "data")
# Use 'celltype2' for labels
reference_labels <- pdac_reference_seurat$celltype2

# Build your SingleCellExperiment object with detailed labels
reference_sce <- SingleCellExperiment(
  assays = list(logcounts = reference_expression),
  colData = DataFrame(label.main = reference_labels)
)

# Prepare dataset
# Run if it did not run before
#merged_seurat.harmony <- JoinLayers(merged_seurat.harmony)
query_expression <- GetAssayData(merged_seurat.harmony, assay = "RNA",layer = "data")
query_sce <- SingleCellExperiment(
  assays = list(logcounts = query_expression),
  colData = merged_seurat.harmony@meta.data
)

# Run SingleR annotation
singleR_annotations <- SingleR(
  test = query_sce,
  ref = reference_sce,
  labels = reference_sce$label.main
)


# Save the Seurat object after doublet detection
#saveRDS(singleR_annotations, file = "singleR_annotations.rds")
# Load the Seurat object after doublet detection
#singleR_annotations <- readRDS("singleR_annotations.rds")

#Add SingleR labels to Seurat object
merged_seurat.harmony$SingleR_Labels <- singleR_annotations$labels

#UMAP
DimPlot(merged_seurat.harmony, reduction = "umap", group.by = "SingleR_Labels", label = TRUE, repel = TRUE, pt.size = 0.5)
DimPlot(merged_seurat.harmony, reduction = "umap", split.by = "group", group.by = "SingleR_Labels", label = TRUE, repel = TRUE, pt.size = 0.5)

# 11. Identify DEGs between tumor and normal pancreatic ductal cells   ---------

# Set cell identities to SingleR-predicted cell types
Idents(merged_seurat.harmony) <- "SingleR_Labels"

# Create a new group combining label and condition
merged_seurat.harmony$Cell_Group <- paste(
  merged_seurat.harmony$SingleR_Labels,
  merged_seurat.harmony$group,
  sep = "_"
)

# Set identity to the combined cell type and condition label
Idents(merged_seurat.harmony) <- "Cell_Group"

# Differential expression: tumor ductal vs normal ductal
ductal_DEGs <- FindMarkers(
  merged_seurat.harmony,
  ident.1 = "pancreatic ductal cell_tumor",
  ident.2 = "pancreatic ductal cell_normal",
  logfc.threshold = 1,
  min.pct = 0.5,
  test.use = "MAST",
)

# Save the DEGs to a CSV file
write.csv(ductal_DEGs, file = "ductal_DEGs.csv")

# RNF149 division
# 1. Subset to ductal tumor cells
ductal_tumor_cells <- subset(seurat_obj, idents = "pancreatic ductal cell_tumor")

# 2. Extract raw counts of RNF149 in ductal tumor cells
RNF149_counts <- GetAssayData(ductal_tumor_cells, slot = "counts")["RNF149", ]

# 3. Identify cells with RNF149 expression > 0
nonzero_cells <- names(RNF149_counts[RNF149_counts > 0])

# 4. Compute median expression (non-zero only)
median_RNF149 <- median(RNF149_counts[nonzero_cells])

# 5. Create a group label: High, Low, Zero
ductal_tumor_cells$RNF149_group <- "Zero"  # default
ductal_tumor_cells$RNF149_group[nonzero_cells] <- ifelse(
  RNF149_counts[nonzero_cells] >= median_RNF149,
  "RNF149 High",
  "RNF149 Low"
)

# 6. Set identity to the RNF149 group
Idents(ductal_tumor_cells) <- ductal_tumor_cells$RNF149_group

# 7. Check group sizes
table(Idents(ductal_tumor_cells))

FeaturePlot(merged_seurat.harmony, features = c("RNF149"), split.by = "group")

VlnPlot(
  merged_seurat.harmony,
  features = "RNF149",
  split.by = "group",
  pt.size = 0
) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("RNF149 Expression Across Cell Types")


# MBOAT7 division
# Subset the tumor pancreatic ductal cells
tumor_ductal <- subset(merged_seurat.harmony, idents = "pancreatic ductal cell_tumor")

# 1. Extract raw counts for MBOAT7
MBOAT7_counts_all <- GetAssayData(tumor_ductal, slot = "counts")["MBOAT7", ]

# 2. Identify cells with expression > 0 (non-zero)
cells_with_MBOAT7 <- names(MBOAT7_counts_all[MBOAT7_counts_all > 0])

# 3. Calculate median of NON-ZERO cells
median_MBOAT7 <- median(MBOAT7_counts_all[cells_with_MBOAT7])

# 4. Classify NON-ZERO cells into High/Low based on median
tumor_ductal$MBOAT7_group <- "Zero"  # Default group
tumor_ductal$MBOAT7_group[cells_with_MBOAT7] <- ifelse(
  MBOAT7_counts_all[cells_with_MBOAT7] >= median_MBOAT7,
  "MBOAT7 High",
  "MBOAT7 Low"
)

# 5. Verify group distribution
table(tumor_ductal$MBOAT7_group)

# 6. Set identities
Idents(tumor_ductal) <- tumor_ductal$MBOAT7_group

FeaturePlot(merged_seurat.harmony, features = c("MBOAT7"), split.by = "group")

VlnPlot(
  merged_seurat.harmony,
  features = "MBOAT7",
  split.by = "group",
  pt.size = 0
) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("MBOAT7 Expression Across Cell Types")
