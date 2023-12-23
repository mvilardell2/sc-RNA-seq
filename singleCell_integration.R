###############################################################################
########## integrate scRNA-Seq datasets to correct for batch effects ##########

setwd("C:/Users/MARINA/Documents/BIOINFORMATICS/scRNA_seq")

# ------- load libraries ------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# ------ get data location -------
# We need the location of the datasets that we want to integrate
dirs <- list.dirs(path = 'datasets/', recursive = F, full.names = F) # 7 folders


# ------ Create the Seurat object -------
# loop to create the Seurat objects with each of the matrices for each dataset --- in total 7 Seurat Objects.
for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  cts <- ReadMtx(mtx = paste0('datasets/',x,'/matrix.mtx.gz'),
          features = paste0('datasets/',x,'/features.tsv.gz'),
          cells = paste0('datasets/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}

# To analize each of the datasets separately would be time consuming, therefore we will merge them. 


# ------- merge datasets ---------
merged_seurat <- merge(HB17_background, y = c(HB17_PDX, HB17_tumor, HB30_PDX, HB30_tumor, HB53_background,
                             HB53_tumor),
      add.cell.ids = ls()[3:9],
      project = 'HB')

merged_seurat

 
# --------- QC & filtering ------------

View(merged_seurat@meta.data)

# create a sample column that tells the patients and the tissue from the barcode
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column into Patient, type and barcode (separator='_')
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Patient', 'Type', 'Barcode'), 
         sep = '_')

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mitoPercent < 10)

# Number of cells before filtering: 
merged_seurat

# Number of cells after filtering:
merged_seurat_filtered



# ------------- perform standard workflow steps to figure out if we see any batch effects --------
# Normalization
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)

# Find features
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)

# Scale data
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)

# Run PCA
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered) #Dimension of the dataset
# In this case, we utilize all 20 dimensions

# Find Neighbors
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)

# Clustering
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)


merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


# ------------ Visualization ------------------
# Color the cells by patient
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient')

# Color the cells by tissue type
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type',
        cols = c('red','green','blue'))

grid.arrange(p1, p2, ncol = 2, nrow = 2)
# Differences may be due to technical effects rather than biological, so we need to correct for batch effects


# ------------------ perform integration to correct for batch effects ----------------
# Split the Seurat object by Patient, to see the batch effects
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Patient') #We have a List with the objects for each patient
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                       anchor.features = features, reduction = 'rpca')

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)
# Now we can run a single integrated analysis on all cells!

# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)


p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)








