Name: Meiheng Liang
Programming Language: [R]
Date: [Mat 2023]
Description:
This script is a demonstration the application of Cell cell communication analysis using scRNA data from patients suffered from atherosclerotic plaque. Seurat pipeline and CellChat package was used in analyzing cellular composition of tissue within calcified core and proximal region. The data which available from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677

############################################################################
Required files:
Patient 1 AC: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4837523
Patient 1 PA: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4837524

############################################################################
Required packages and library:
##load packages
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages("devtools")
devtools::install_github("sqjin/CellChat")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("Seurat")
##Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(rhdf5)
library(ggplot2)
library(CellChat)

############################################################################

Following instruction within the script for code execution.

For cellchat, detailed documentation are available:
https://github.com/jinworks/CellChat/tree/main/tutorial

############################################################################
# Output files: 
ac_top30_gene.csv
pa_top30_gene.csv
