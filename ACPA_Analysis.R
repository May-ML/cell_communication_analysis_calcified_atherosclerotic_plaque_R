
# package prep ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("sqjin/CellChat")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("Seurat")
install.packages("ggplot2")
install.packages("remotes")
remotes::install_github("sqjin/CellChat")
tools::package_dependencies("ggplot2", recursive = TRUE, which = "Imports")
##Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(rhdf5)
library(ggplot2)
library(CellChat)

# load data ---------------------------------------------------------------

#Check if the file is in your working directory
dir()
#Change your working directory to location of files
setwd("C:/Users/Chris/git-test/528Project")
#Load in 10x Genomics data
AC.data <- Read10X("C:/path/to/your/GSM4837523_02dat20190515tisCARconDIS_featurebcmatrixfiltered/outs/filtered_feature_bc_matrix")
PA.data <- Read10X("C:/path/to/your/GSM4837524_01dat20190515tisCARconHEA_featurebcmatrixfiltered/outs/filtered_feature_bc_matrix")
#Create RNA object, where the counts
ac <- CreateSeuratObject(counts = AC.data, project = "atherosclerotic plaque ", min.cells = 0, min.features = 200)
#Create RNA object, where the counts
pa <- CreateSeuratObject(counts = PA.data, project = "Proximal Region control ", min.cells = 0, min.features = )
#standardized min.cell=3, and min.feature=200, reference used min.cell=0, min.feature=0 
ac[["percent.MT"]] <- PercentageFeatureSet(ac, pattern = "^MT-")
pa[["percent.MT"]] <- PercentageFeatureSet(pa, pattern = "^MT-")

#Carotid artery atherosclerotic core # Features, samples, 1 assay
ac #check data
pa


# Data preprocessing AC ---------------------------------------------------

head(ac) #check data structure
# Use character vector names directly
VlnPlot(ac, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
# Display these column types 
ac1 <- subset(ac, subset = nFeature_RNA > 200 & nCount_RNA < 4000 & percent.MT < 10)

VlnPlot(ac1, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

##normalization
#pre_norm
# Highly variable gene shown pre-normalization
VlnPlot(ac1,features="ANXA1") #immune cell marker, or ACTA2 a marker for smooth muscle cell
ac1 <- NormalizeData(ac1, normalization.method = "LogNormalize", scale.factor = 10000) # Scale factor is the magnitude of your normalization. scale factor referred to original document 
# Highly variable gene shown pre-normalization
VlnPlot(ac1,features="ANXA1") # Look at the difference after normalization
acmedgene<-round(median(ac1$nFeature_RNA))
## Identification key variable 
# adjust thereshold as necessary
ac1 <- FindVariableFeatures(ac1, selection.method = "vst", nfeatures = 979)  # Identifies features that are outliers on a 'mean variability plot'.

## Scaling the data  #scaling for 10,000 here
all.genes <- rownames(ac1)
ac1 <- ScaleData(ac1, features = all.genes)

# Linear Dimension Reduction ----------------------------------------------

# principal component 1,2,3 decenting importance in impacting data
ac1 <- RunPCA(ac1, features = VariableFeatures(object = ac1))
VizDimLoadings(ac1, dims = 1:2, reduction = "pca")
DimPlot(ac1, reduction = "pca")

DimHeatmap(ac1, dims = 1, cells = 500, balanced = TRUE)
##principal component determines good separation between upregulated or downregulated genes, if the separation of gene expression in heatmap is not clear, wrong with principal component
##ideally, heat map contain top expression gene and least expression gene in control versus disease, determine a separation point/cutoff point for PCA based on graph visualization

DimHeatmap(ac1, dims = 1:5, cells = 500, balanced = TRUE) #clearity dropped after 5 at first DimPlot
DimHeatmap(ac1, dims = 1:15, cells = 500, balanced = TRUE) 


# Determine Dimensionality ------------------------------------------------
ac1 <- JackStraw(ac1, num.replicate = 100)
ac1 <- ScoreJackStraw(ac1, dims = 1:20)
ElbowPlot(ac1,ndims = 30)   ##looking for highly variated genes, when the elbow tilting toward 0 means less variation and less significant

## Cells Clustering
# Use elbow plot to determine how many dimensions to look at, 
# Here we can go to 1:25 but even 1:40 works
ac2 <- FindNeighbors(ac1, dims = 1:25) #
ac2 <- FindClusters(ac2, resolution = 1.0) #resolution adjust according to cell group, higher the resolution smaller number of cell in scope, more detail for each
# Think of resolutions like microscopes. You can use different resolutions
# to change how "zoomed in" you want to be



# Non-linear dimensional reduction (UMAP/tSNE) ----------------------------

ac2@active.ident # Looking at the current amount of clusters 
ac2 <- RunUMAP(ac2, dims = 1:16)  # of cluster from previous run
DimPlot(ac2, reduction = "umap", label = TRUE)
# Looking at one specific gene
FeaturePlot(ac2, features = c("NKG7")) 

# saveRDS(ac2,file = "AC Clustering.rds")

#visualize expression profile
FeaturePlot(ac2, features = c("TAGLN")) 


# Identify differentially expressed features and Annotation ---------------
# Get all genes that are positive + assign markers
ac_gene <- FindAllMarkers(ac2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Count the number of cells in each cluster
ac_cluster_size <- table(Idents(ac2))
print(ac_cluster_size)

# Ordering the markers by logFC and getting the top 20 for each cluster
ac_top_markers <- ac_gene %>% 
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)  

write.csv(ac_top_markers, "ac_top30_gene.csv")
# Viewing the top markers
head(ac_top_markers)
# Then go to https://panglaodb.se/samples.html and manually look up genes, 
# Look at each cluster, significant genes (p-value), 
# Rename clusters ID  - step for annotation
new.cluster.ids <- c("T cell", "Cytotoxic T cells","Cytotoxic T cells", "Macrophages", "Cytotoxic T cells", "NK cell", 
                     "NK cell", "basal cell", "Cytotoxic T cells", "Endothelial cells", "Smooth muscle cells",
                     "Macrophages", "Neutrophils", "Mast cells", "endothelial cells", "B cells", "monocytes")
levels(ac1)
names(new.cluster.ids) <- levels(Idents(ac2))  # Use Idents to get the levels of the current identities
ac2 <- RenameIdents(ac2, new.cluster.ids)

#visualize via vln plot: 

# Plot the Umap but with the annotated clusters
DimPlot(ac2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


VlnPlot(ac2, features = c("SPP1", "SFRP5", "IBSP", "CRTAC1","ITLN1","DKK2", "F5", "FN1"), ncol = 8)


# Proximal Region ---------------------------------------------------------

# Preprocessing -----------------------------------------------------------

# Use character vector names directly
VlnPlot(pa, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

# Display these column types 
pa1 <- subset(pa, subset = nFeature_RNA > 200 & nCount_RNA < 4000 & percent.MT < 10)

VlnPlot(pa1, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

# Highly variable gene shown pre-normalization
VlnPlot(pa1,features="HLA-DPB1") #immune cell marker, or ACTA2 a marker for smooth muscle cell


# Normalization -----------------------------------------------------------
pa2 <- NormalizeData(pa1, normalization.method = "LogNormalize", scale.factor = 10000) # Scale factor is the magnitude of your normalization. scale factor referred to original document 
# Highly variable gene shown pre-normalization
VlnPlot(pa2,features="HLA-DPB1") # Look at the difference after normalization
pamedgene<-round(median(pa2$nFeature_RNA))
## Identification key variable 
# adjust thereshold as necessary
pa2 <- FindVariableFeatures(pa2, selection.method = "vst", nfeatures = 1186)  # Identifies features that are outliers on a 'mean variability plot'.

## Scaling the data  #scaling for 10,000 here
all.genes <- rownames(pa2)
pa2 <- ScaleData(pa2, features = all.genes)


# linear dimensional reduction --------------------------------------------

# principal component 1,2,3 decensing importance in impacting data
pa2 <- RunPCA(pa2, features = VariableFeatures(object = pa2))
VizDimLoadings(pa2, dims = 1:2, reduction = "pca")
DimPlot(pa2, reduction = "pca")

DimHeatmap(pa2, dims = 1, cells = 500, balanced = TRUE)
##principal component determines good separation between upregulated or downregulated genes, if the separation of gene expression in heatmap is not clear, wrong with principal component
##ideally, heat map contain top expression gene and least expression gene in control versus disease, determine a separation point/cutoff point for PCA based on graph visualization

DimHeatmap(pa2, dims = 1:2, cells = 500, balanced = TRUE)
DimHeatmap(pa2, dims = 1:15, cells = 500, balanced = TRUE) 


# Determine dimensionality -----------------------------------------------
pa2 <- JackStraw(pa2, num.replicate = 100)
pa2 <- ScoreJackStraw(pa2, dims = 1:20)
ElbowPlot(pa2,ndims = 50)   ##looking for highly variated genes, when the elbow tilting toward 0 means less variation and less significant

## Cells Clustering
# Use elbow plot to determine how many dimensions to look at, 
pa2 <- FindNeighbors(pa2, dims = 1:7) #
pa2 <- FindClusters(pa2, resolution = 1.0) #resolution adjust according to cell group, higher the resolution smaller number of cell in scope, more detail for each
# Think of resolutions like microscopes. You can use different resolutions
# to change how "zoomed in" you want to be


# Visualization --------------------------------------------------------------------
pa2@active.ident # Looking at the current amount of clusters 
pa2 <- RunUMAP(pa2, dims = 1:8)  # of cluster from previous run
DimPlot(pa2, reduction = "umap", label = TRUE)
# Looking at one specific gene
FeaturePlot(pa2, features = c("HLA-DPB1")) 


# Annotation --------------------------------------------------------------
# Get all genes that are positive + assign markers
pa_gene <- FindAllMarkers(pa2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(pa_gene)
# Count the number of cells in each cluster
pa_cluster_size <- table(Idents(pa2))
print(pa_cluster_size)

# Ordering the markers by logFC and getting the top 20 for each cluster
pa_top_markers <- pa_gene %>% 
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)  

write.csv(pa_top_markers, "pa_top30_gene.csv")
# Viewing the top markers

head(pa_top_markers)
#levels(pa2)

new.cluster.ids <- c("Endothelial cells", "Smooth muscle cells","Fibroblasts", "T cells", "Cytotoxic T cells", "Endothelial cells", 
                     "macrophages", "metastatic epithelial cell", "Neutrophils")

names(new.cluster.ids) <- levels(Idents(pa2))  # Use Idents to get the levels of the current identities
pa2 <- RenameIdents(pa2, new.cluster.ids)

#visualize via vln plot: 
VlnPlot(pa2, features = c("APOD", "MFAP5", "PLA2G2A", "C3"), ncol = 4)
VlnPlot(pa2, features = c("IL6", "MLPH", "HLA-DQA1", "ACKR1"), ncol = 4)

# Plot the Umap but with the annotated clusters
DimPlot(pa2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# cc communication --------------------------------------------------------

# Print a summary of the Seurat object
pa2

# Check the first few rows of meta.data to confirm the update
head(pa2$meta.data)

metadata <- pa2@meta.data
# Suppose you have 9 clusters and you've determined each cluster's cell type
cluster_names <- c("Endothelial cells", "Smooth muscle cells","Fibroblasts", "T cells", "Cytotoxic T cells", "Endothelial cells", 
                   "macrophages", "metastatic epithelial cell", "Neutrophils")

# Assume cluster numbers in Seurat are 0 to 8 (adjust if your clusters start from 1)
# Mapping cluster numbers to names
names(cluster_names) <- 0:8  # Adjust this range based on your actual cluster numbers

# Update the metadata to include cell type names
metadata$cell_type <- cluster_names[as.character(pa2@meta.data$seurat_clusters)]

# Check the first few rows to ensure the update
head(metadata)


# Create the CellChat object
cellchat <- createCellChat(object = expression_matrix, meta = metadata, group.by = "seurat_clusters")

# Load the default signaling database, adjust as per organism (human, mouse, etc.)
data(CellChatDB.human)
cellchat@DB <- CellChatDB.human

# Proceed with your analysis, e.g., identify and analyze cell-cell communication networks
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat)
cellchat <- computeCommunProb(cellchat)

# Analysis and visualization of networks
netAnalysis(cellchat)
netVisual(cellchat)





























