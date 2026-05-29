# Meningeal-lymphatics-as-a-gateway-for-extracellular-vesicle-brain-delivery

This repository contains the standalone R scripts used for the analysis in the manuscript "Meningeal-lymphatics-as-a-gateway-for-extracellular-vesicle-brain-delivery".

## 1. System Requirements

### Hardware requirements
This script requires a standard computer with enough RAM to support the in-memory operations. 
* **Minimum:** 32 GB RAM
* **Recommended for full dataset:** 125+ GB RAM (or HPC environment)

The scripts have been specifically tested on the following system configuration:
* **OS:** Linux Rocky Linux release 8.10 (Green Obsidian)
* **RAM:** 128 GB
* **CPU:** AMD Ryzen Threadripper PRO 5995WX 64-Cores

### Software requirements
#### R Dependencies
This software has been tested on **R version 4.5.1**.
The core scripts require the following R packages. The versions the software has been tested on are listed below:

* `readxl` (v1.4.5)
* `dplyr` (v1.1.4)
* `tidyr` (v1.3.1)
* `clusterProfiler` (v4.16.0)
* `org.Hs.eg.db` (v3.21.0)
* `ggplot2` (v4.0.1)
* `circlize` (v0.4.16)
* `tidyverse` (v2.0.0)
* `ggalluvial` (v0.12.5)

## 2. Installation Guide

### Instructions
The required R packages can be installed via CRAN and Bioconductor. Open your R console and run the following commands to set up the environment.

```R
# 1. Install CRAN packages
cran_packages <- c("readxl", "dplyr", "tidyr", "tidyverse", "ggplot2", "ggalluvial", "circlize")
install.packages(cran_packages)

# 2. Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

bioc_packages <- c("clusterProfiler", "org.Hs.eg.db")
BiocManager::install(bioc_packages)


## 3 & 4. Instructions for Use and Reproduction
1. Data Source (For Reviewers & Users)
Raw LC/MS Data (PRIDE):
The raw LC-MS/MS proteomics data associated with this manuscript have been deposited to the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD079034.

Note to Reviewers: The PRIDE dataset is currently on hold until publication. To access the raw data during the peer review process, please log in to the PRIDE website using the following details:

Project accession: PXD079034

Token: yW9qu0IJw9iu

The complete raw proteomic expression dataset containing the HEK and MSC samples is directly provided within this GitHub repository. Reviewers and users can utilize EV proteomic Hek vs MSC.xlsx as the primary raw data file to fully replicate and run the analytical pipeline.

2. Instructions to Run
Clone this repository to your local machine.
Ensure the raw Excel data file (EV proteomic Hek vs MSC.xlsx) is correctly placed in your designated working directory.
Open the main analysis script (MSC-EV proteomic analysis for submission.R) in RStudio.
Under the 1. Data Loading section, modify the following directory path to match where you saved the downloaded dataset:
file_path <- "your directory/EV proteomic Hek vs MSC.xlsx"
Execute the script sequentially from top to bottom.

3. Expected Output
Running the script will process the Common, MSC Unique, and HEK Unique proteins, map the Ligand-Receptor-Pathway axis, and generate the quantitative figures used in the manuscript, including:
Dot Plot: "Functional Profile of MSC Exosomes (Total)" displaying enriched KEGG pathways sorted by statistical significance.
Chord Diagram: Visualizing enriched pathways (e.g., Endocytosis, PI3K-Akt, Regulation of actin cytoskeleton) and target genes.
Sankey/Alluvial Diagram: "Mechanism of MSC Exosome Uptake (Balanced View)" linking Ligand Families to Lymphatic Receptors and Intracellular Pathways.
CSV Output: MSC_Exosome_Ligand_Candidates_87.csv containing the identified ligand candidate list.

4. Expected Run Time
Processing the datasets and running downstream enrichment tools requires substantial computational resources.
Expected run time: Approximately 30 minutes on a standard machine with at least 32GB RAM.
