##### script to perform standard workflow steps to analyze single cell RNA-Seq data #####
#########################################################################################

# data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1
# data source: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0         

setwd("C:/Users/MARINA/Documents/BIOINFORMATICS/scRNA_seq")

# load libraries
library(Seurat)
library(tidyverse)

# Load the NSCLC dataset
nsclc.sparse.m <- Read10X_h5(filename = 'Datasets/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
str(nsclc.sparse.m) 
cts <-  nsclc.sparse.m$`Gene Expression` #Not a seurat object yet
cts[1:10,1:10] #In the rows we have the features(genes), and in the columns the cell barcode


# Initialize the Seurat object with the raw (non-normalized data).
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200) 
#min.cells->keep all the features that have at least expression in 3 cells
#min.features->keep the cells that have at least 200 features
str(nsclc.seurat.obj)
nsclc.seurat.obj
# 29552 features across 42081 samples


# 1. --------- QC -----------
# QC metrics are found in @meta.data slot

# QC metrics: 
# - Number of unique genes detected in each cell
# - Total number of molecules detected within a cell
# - % of reads that map to the mitochondrial genome

View(nsclc.seurat.obj@meta.data)

# % MT reads-->Calculate the mitochondrial genes. Low quality cells will have high expression
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

# Violin plot
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Lot of cells have a high percentage of mitochondrial reads, and these need to be removed (cells of not good quality)

FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm') #A good quality dataset should follow the straight line. 



# 2.--------- Filtering low quality cells-----------------

# We filter out those cells that: 

# - Have feature counts over 2500 and less than 200
# - Have more than 5% f mitochondrial counts

nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)
# 29552 features across 24708 samples


# 3.---------- Normalize data ----------
#In order to compare between different cells
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)


# 4.---------- Identify highly variable features --------------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


# 5.------ Scaling -------------
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)


# SLOTS: 
# - Counts: raw data
# - Data: logs normalized counts 
# - Scale.data: after the scale


# 6.------------ Perform Linear dimensionality reduction (PCA) --------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))
# We can see the genes that have a positive and negative score for the PC_1, 2, 3 ...

# visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)

#Heatmap
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)
# These are the genes that exhibit heterogeneity

# determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj)
#Consider the PC that explain a high percentage of variability after which the variance of the PC do not very much. Better to consider more PC than lower number of PC.


# 7.-------- Clustering ------------
#Cluster cells that have similar expression pattern
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15) #Dims parameter is the number of PC that we have selected previously.

# understanding resolution -  assign cells to clusters
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))

# For each resolution, there is a new column created in the metadata, that tells how many clusters are there for each resolution.
View(nsclc.seurat.obj@meta.data)

# Which is the best resolution?
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
# With resolution 0.1 we identify 8 clusters, and seems that there are separated.
# With resolution 0.5 there are 12 clusters, and some of them are overlaped.
# In this case, we choose resolution 0.1. But if we have a larger dataset, we may require larger resolution to separate the cells

# setting identity of clusters
head(Idents(nsclc.seurat.obj),5)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
head(Idents(nsclc.seurat.obj)) # Now, the cells are identified by 8 clusters



# ------------------ non-linear dimensionality reduction --------------
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label

# individual clusters
DimPlot(nsclc.seurat.obj, reduction = "umap")
# We have 8 clusters (resolution 0.1)




