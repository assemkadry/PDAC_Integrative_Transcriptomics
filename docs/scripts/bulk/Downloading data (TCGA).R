# Install and load required packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

# Define query parameters
query <- GDCquery(
  project = "TCGA-PAAD",  # Pancreatic adenocarcinoma
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",  # STAR counts
  platform = "Illumina",  # Illumina platform
  access = "open",  # Open-access data
  sample.type = c("Primary Tumor", "Metastatic")  # Include both primary and metastatic samples
)

# Download data
GDCdownload(query)

# Prepare and load data
rna_data <- GDCprepare(query)

library(SummarizedExperiment)

# View first few rows of the expression matrix
head(assay(rna_data))

# Save as CSV file
saveRDS(rna_data, "TCGA_PAAD_RNAseq_STAR_Counts.rds")

TCGA_PAAD_data <- readRDS("TCGA_PAAD_RNAseq_STAR_Counts.rds")

library(dplyr)
library(here)

colData(TCGA_PAAD_data) %>%
  colnames() %>%
  tail()
