# 🔬 Integrated Analysis of Bulk and Single-Cell Transcriptomics in Pancreatic Cancer

This repository accompanies the study:

**_Integrative Analysis of Bulk and Single-Cell Transcriptomics Reveals Two Novel Genes Associated with Poor Prognosis in Pancreatic Cancer_**

> 🧪 This study is currently under review. Results and interpretations will be available upon publication.

---

## 📚 Table of Contents

- [Overview](#-overview)
- [Workflow Diagram](#-workflow-diagram)
- [Repository Structure](#-repository-structure)
- [Data Sources](#-data-sources)
- [Setup](#️-setup-and-requirements)
- [Outputs](#-outputs)
- [Author](#-author)
- [Keywords](#-keywords)
- [License](#-license)

---

## 🧪 Overview

This project integrates single-cell and bulk RNA-seq datasets to identify molecular drivers of pancreatic ductal adenocarcinoma (PDAC). Two genes, RNF149 and MBOAT7, were identified through comparative expression in ductal cells and analyzed further in bulk TCGA data to understand their role in immune and lipid metabolic remodeling in PDAC.

---

## 🔬 Flow of Analysis

1. **Single-cell RNA-seq (GSE155698)**: Clustering, ductal cell annotation, DEA between tumor and normal ductal cells.
2. **Gene selection**: RNF149 (upregulated) and MBOAT7 (downregulated) from single-cell analysis.
3. **Bulk RNA-seq (TCGA-PAAD)**: DEA based on RNF149/MBOAT7 expression levels.
4. **Pathway & interactome analysis**: KEGG, Reactome, GO, GSEA, and protein interaction analysis.

---

## 🧭 Workflow Diagram

📌 _Coming Soon: Full visual pipeline diagram._

---

## 📁 Repository Structure

```
docs/
├── scripts/
│   ├── scRNA/
│   │   └── scRNA analysis.R
│   ├── bulk/
│   │   ├── Downloading data (TCGA).R
│   │   ├── Preprocessing + DEA (RNF149).R
│   │   └── Preprocessing + DEA (MBOAT7).R
│   └── enrichment/
│       ├── RNF149 PATHWAY ENRICHMENT (KEGG + REACTOME ).R
│       ├── MBOAT7 Pathway Enrichment (Kegg + Reactome).R
│       ├── RNF149-INTACT-Multi- INTERACTORS-DEGS-GO-ENRICHMENT.R
│       └── Ranked_GSEA.r
```

---

## 📂 Data Sources

- 🔸 **Single-cell RNA-seq**: GSE155698 (PDAC and adjacent normal tissue)  
- 🔸 **Bulk RNA-seq**: TCGA-PAAD (via TCGAbiolinks)

---

## ⚙️ Setup and Requirements

Install required R packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# CRAN packages
install.packages(c(
  "tidyverse", "data.table", "ggplot2", "pheatmap", "cowplot"
))

# Bioconductor packages
BiocManager::install(c(
  "Seurat", "DESeq2", "org.Hs.eg.db", "AnnotationDbi", 
  "clusterProfiler", "ReactomePA", "DOSE", "msigdbr", 
  "TCGAbiolinks", "SummarizedExperiment", "fgsea"
))
```

---

## 📊 Outputs

- 🔬 **Single-cell**: UMAP plots, DEGs between ductal tumor vs. normal cells.
- 🧬 **Bulk RNA-seq**: DEGs based on RNF149/MBOAT7 expression, volcano plots.
- 🧠 **Enrichment**: GO/KEGG/Reactome results for significant genes.
- 📈 **GSEA**: Ranked enrichment analysis.

---

## 👨‍💻 Authors

**Assem K. Elsherif**
**Nourine Mamdouh Sabry Abdelfattah**
**Sajda Hussien Salah Tahoun**
**Sondos Ameen El-Sayed Mohammed Awad**
**Moaz Mohamed ElShiekh**

---

## 🧠 Keywords

Pancreatic Cancer · Bulk RNA-seq · Single-cell RNA-seq · RNF149 · MBOAT7 · TCGA · Seurat · DESeq2 · Enrichment

---

## 📜 License

This project is licensed under the MIT License.
