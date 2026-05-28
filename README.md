# Meningeal-lymphatics-as-a-gateway-for-extracellular-vesicle-brain-delivery
Raw Mouse Meninges data and analysis code
# [Project Title: e.g., Mouse Hippocampus Aging scRNA-seq Analysis]

This repository contains the standalone R scripts and a demo dataset used for the analysis in the manuscript "[Your Manuscript Title]".

## 1. System Requirements

### Hardware requirements
This script requires a standard computer with enough RAM to support the in-memory operations. 
* **Minimum:** 32 GB RAM
* **Recommended for full dataset:** 125+ GB RAM (or HPC environment)

### Software requirements
The scripts have been tested on the following operating systems:
* Linux: [Rocky Linux release 8.10 (Green Obsidian)]
* RAM: 128 GB
* CPU: AMD Ryzen Threadripper PRO 5995WX 64-Cores

**Tested R version:** R version [v4.5.1]

#### R Dependencies
This software has been tested on **R version 4.5.1**.
The core scripts require the following R packages. The versions the software has been tested on are listed below:

* `Seurat` (v5.4.0)
* `dplyr` (v1.1.4)
* `ggplot2` (v4.0.1)
* `patchwork` (v1.3.2)
* `harmony` (v1.2.4)
* `celldex` (v1.20.0)
* `SingleR` (v2.12.0)
* `scales` (v1.4.0)
* `ggrepel` (v0.9.6)
* `org.Mm.eg.db` (v3.22.0)
* `stringr` (v1.6.0)
* `clusterProfiler` (v4.18.4)
* `GOplot` (v1.0.2)
* `tidyr` (v1.3.2)
* `enrichplot` (v1.30.4)
* `circlize` (v0.4.17)
* `ComplexHeatmap` (v2.26.0)
* `RColorBrewer` (v1.1-3)

## 2. Installation Guide

### Instructions
The required R packages can be installed via CRAN and Bioconductor. Open your R console and run the following commands to set up the environment.

```R
# 1. Install CRAN packages
cran_packages <- c("Seurat", "dplyr", "ggplot2", "patchwork", "harmony","scales", "ggrepel", "stringr", "tidyr", "circlize", "RColorBrewer", "GOplot")
install.packages(cran_packages)

# 2. Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")}

bioc_packages <- c("celldex", "SingleR", "org.Mm.eg.db", "clusterProfiler", "enrichplot", "ComplexHeatmap")
BiocManager::install(bioc_packages)

## 3 & 4. Instructions for Use and Reproduction

Due to the large size of the raw datasets, we do not provide a subsetted demo dataset in this repository. Instead, users and reviewers can fully reproduce the manuscript's results by downloading the complete datasets from public repositories.

### 1. Data Download (For Reviewers)
The single-cell RNA sequencing data for the mouse meninges (Young, Old, and Old-EV groups) has been deposited in the Gene Expression Omnibus (GEO) under the accession number **[GSE333240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333240)**.

> **Note to Reviewers:** The GEO dataset is currently on hold until publication. To access the raw data, please visit the link above and enter the private access token: **`szcbmaiohtobvwb`**

Please download the `barcodes.tsv`, `features.tsv`, and `matrix.mtx` files for the following samples and place them in the `data/raw/` directory:
* GSM9759670: Young Mouse Meninges
* GSM9759671: Old Mouse Meninges
* GSM9759672: Old-EV Mouse Meninges

*(Note: LC/MS data from PRIDE can be downloaded and placed in the `data/raw/` directory following similar steps. Accession: [Insert PRIDE Accession])*

### 2. Instructions to Run
1. Clone this repository to your local machine or HPC environment.
2. Ensure all downloaded GEO and PRIDE data files are correctly placed in the `data/raw/` folder.
3. Open the main analysis script (e.g., `main_analysis.R`) in RStudio.
4. Execute the script sequentially. 

### 3. Expected Output
Running the script will generate the following outputs in the `results/` folder:
* Quality control plots
* UMAP dimension reduction visualizations (e.g., cell type clustering including Endothelial cells)
* Differentially Expressed Genes (DEG) tables and Gene Ontology (GO) pathway enrichment results.

### 4. Expected Run Time
Processing the full dataset requires substantial computational resources.
* **Expected run time:** Approximately [e.g., 1-2 hours] on a standard machine with 32GB+ RAM.
