##############  Analysis and visualiation of spatial datasets with Seurat #####################
###############################################################################################

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


setwd("C:/Users/MARINA/Documents/BIOINFORMATICS/scRNA_seq")

### We will analyse a dataset generated with Visium 10x Genomics. 
# This data sets is from mous brain slices

InstallData("stxBrain")
?stxBrain

brain <- LoadData("stxBrain", type = "anterior1")
str(brain)

# -------- Data preprocessing -----------
# Explore counts
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
# These plots demonstrate the variance in the molecular count across the spots, that will depend on the tissue anatomy.

## NORMALIZATION 
# Therefore, in order to normalize the data, we will use SCTransform, which detects high-variance features, and sotres the data in the SCT assay. 
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)


# ------ Gene expression visualization ----------
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))


# --------- Dimension reduction and clustering -------------
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2





