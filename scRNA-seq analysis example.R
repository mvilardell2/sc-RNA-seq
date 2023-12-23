################################################################################
#################   SINGLE CELL RNA-SEQ ANALYSIS EXAMPLE   #####################
################################################################################


#data: Peripheral Blood Mononuclear Cells (PBMC) from 10X Genomics. 

setwd("C:/Users/MARINA/Documents/BIOINFORMATICS/scRNA_seq")

# Load libraries
library(dplyr)
library(Seurat)

# ----------- Data Preparation --------------
# Load the dataset
pbmc.data<-Read10X(data.dir= "Datasets/filtered_gene_bc_matrices/hg19/") #Not a Seurat Object
str(pbmc.data)
pbmc.data[1:10,1:10]


# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc


# What does data in a count matrix look like?
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]


# -------------- Quality Control-----------------
str(pbmc)
# QC metrics are found in @meta.data slot
head(pbmc@meta.data, 5)
View(pbmc@meta.data)

# Calculate % MT reads
# The [[ operator can add columns to object metadata.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# -------  Filtering low quality cells -----------
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)



# --------- Normalization ----------------
pbmc <- NormalizeData(pbmc)

# Visualize normalized values
head(pbmc[["RNA"]]$data)



# -------- Feature Selection ---------
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# ----------- Scaling data -------------
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# results are stored in: 
pbmc[["RNA"]]$scale.data


# ----------- Linear dimension reduction (PCA) ----------------
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()

# Determine the dimensionality of the dataset 
ElbowPlot(pbmc)


# ----------------- Clustering --------------------

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

View(pbmc@meta.data)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


# -------- non-linear dimensional reduction ----------
pbmc <- RunUMAP(pbmc, dims = 1:10)

# individual clusters
DimPlot(pbmc, reduction = "umap")


# ----- Finding differentially expressed features ------
# The aim of this step is to identify the markers that define a cluster 
# By default, it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. 

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)


##### VISUALIZATION

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# Shows expression of genes across clusters

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD79A"))
#Visualizes feature expression



# ---------- Assigning cell type identity to clusters ---------
#In this dataset, we use canonical markers to easily match the clustering to known cell types
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()




