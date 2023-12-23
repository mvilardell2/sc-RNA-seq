################### script to integrate across conditions using HARMONY #########################
#################################################################################################

setwd("C:/Users/MARINA/Documents/BIOINFORMATICS/scRNA_seq")


# set seed for reproducibility
set.seed(1234)

library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)


# For this example, we are using the 'ifnb' dataset that corresponds to IFNB-stimulated and control PBMCs of human. 

# load dataset
ifnb<-LoadData("ifnb")
str(ifnb)


# -------------QC and filtering -------------
ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = '^MT-')
View(ifnb@meta.data)


# filter
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 & 
                          mito.percent < 5)


# ------------ standard workflow steps -------------
ifnb.filtered <- NormalizeData(ifnb.filtered)
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)
ifnb.filtered <- ScaleData(ifnb.filtered)
ifnb.filtered <- RunPCA(ifnb.filtered)
ElbowPlot(ifnb.filtered)
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')

before <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim')

# Until here we have not performed integration
# Cells grouping by condition are clearly separately. Since these cells are from before and after the treatment, there should be cells similar in both conditions
# SO, we have to integrate to overlay the similar cells



# ----------------- Run Harmony -----------
ifnb.harmony <- ifnb.filtered %>%
  RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

ifnb.harmony@reductions

ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony")
ifnb.harmony.embed[1:10,1:10]
# we can see the harmony components (similar to PCA)


# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize 
after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')
# Cells from both conditions overlay well over each other and do not separate. That is what we wanted to do with integration!


before|after



################### Find markers and cluster identification #########################
#####################################################################################

# script to identify cluster identity 
# Finding markers in every cluster
# Finding conserved markers 
# Finding markers DE between conditions

View(ifnb.harmony@meta.data)

# visualize data; grouped by clusters and by condition
clusters <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
condition <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')

condition|clusters

# So, the data is integrated and then separated in clusters. 
# The goal is to annotate the clusters and to identify the cell types that form the clusters. 

# findAll markers ----------------- This function is more appropriate to use when we have data just from one condition/group
# Identify the cells that form each cluster: 
FindAllMarkers(ifnb.harmony,
               logfc.threshold = 0.25,
               min.pct = 0.1, #test those genes that are found in 0.1 frequency of cells
               only.pos = TRUE, #only positive markers
               test.use = 'DESeq2',
               slot = 'counts')


# In this case we will use findConserved markers

# ----------- findConserved markers -------------

# Notes:
# slot depends on the type of the test used, default is data slot that stores normalized data DefaultAssay(ifnb_harmony) <- 'RNA'

DefaultAssay(ifnb.harmony)

# for example, we want to identify the cluster 3, and compare it with all the other clusters (bc I want to get the markers that are upregulated in cluster 3)
markers_cluster3 <- FindConservedMarkers(ifnb.harmony,
                                         ident.1 = 3,
                                         grouping.var = 'stim')

# it separates cells by condition, and for each compares the cluster 3 with all others, so we get the statistics for each condition
head(markers_cluster3)
# For example, gene FCGR3A, is detected in 97% of the cells in cluster 3 compared to 20% of the cells in other clusters combined.



# let's visualize top features
FeaturePlot(ifnb.harmony, features = c('FCGR3A'), min.cutoff = 'q10')
#min.cutoff= if the expression of the gene for each cell is less than the value, it will be in grey. If the expression is greater, the cells will be colored 

# min-cut off explanation:
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))



#------------ rename cluster 3 ident------------
Idents(ifnb.harmony) #Identify is how our cells are identifies. Untill now, with a number
ifnb.harmony <- RenameIdents(ifnb.harmony, `3` = 'CD16 Mono')

DimPlot(ifnb.harmony, reduction = 'umap', label = T)


# ----- Renam all clusters -----
# cells already have annotations provided in the metadata
View(ifnb.harmony@meta.data)


# Settings cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual anotations (sometimes both)
# need to make sure each cell type forms a separate cluster


# setting Idents as Seurat annotations provided 
Idents(ifnb.harmony) <- ifnb.harmony@meta.data$seurat_annotations
Idents(ifnb.harmony)

DimPlot(ifnb.harmony, reduction = 'umap', label = TRUE)


# -----------Identify Differentially expressed genes (findMarkers) between conditions ---------------------
# Create a new column in the metadata that has the cell type and condition
ifnb.harmony$celltype.cnd <- paste0(ifnb.harmony$seurat_annotations,'_', ifnb.harmony$stim)
View(ifnb.harmony@meta.data)

#Rename the identities
Idents(ifnb.harmony) <- ifnb.harmony$celltype.cnd

DimPlot(ifnb.harmony, reduction = 'umap', label = TRUE)
#So, now the cells belong to two conditions, stimulated and control

##### find markers
# We will compare the cells of CD16 Mono of the control vs the stimulated group: 
b.interferon.response <- FindMarkers(ifnb.harmony, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(markers_cluster3)


FeaturePlot(ifnb.harmony, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')





