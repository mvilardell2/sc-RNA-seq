
# Integration of BCR and GEX data

#################################################################################
#################################################################################

############## BLOC 1. ANALYSE GENE EXPRESSION DATA #############

# ------------ LOAD DATA -----------------

library(hdf5r)
library(Seurat)

# Function to load data and create Seurat object
load_and_create_seurat_object <- function(file_path, project_name) {
  data <- Read10X_h5(file_path, use.names = TRUE, unique.features = TRUE)
  
  # Check and extract 'Gene Expression' modality if present
  if (is.list(data) && "Gene Expression" %in% names(data)) {
    rna_data <- data$`Gene Expression`
  } else if (is.list(data)) {
    rna_data <- data[[1]]
  } else {
    rna_data <- data
  }
  
  # Create Seurat object
  seurat_object <- CreateSeuratObject(counts = rna_data, project = project_name, min.cells = 3, min.features = 100)
  
  return(seurat_object)
}


# Load data
file_path <- "sample_filtered_feature_bc_matrix.h5"
project_name <- "VDJ_GEX_analysis"  
seurat_obj <- load_and_create_seurat_object(file_path, project_name)


# --------- Quality control ------------
head(seurat_obj@meta.data, 5)
# Compute percent mito ratio
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")


VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)


# -------- Filtering low quality cells -----------
table(seurat_obj$orig.ident)
filtered_seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mitoRatio < 10)
table(filtered_seurat_obj$orig.ident)


# ---------- Normalize data ----------
filtered_seurat <- NormalizeData(filtered_seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
head(filtered_seurat[["RNA"]]$data)


# -------- Feature Selection ---------
filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)



# ---------- Scale data -------------
filtered_seurat <- ScaleData(filtered_seurat)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(filtered_seurat), 10)
# Plot variable features with and without labels
var_features <- VariableFeaturePlot(filtered_seurat)
lab_points <- LabelPoints(plot = var_features, points = top10, repel = TRUE)


# ---------- PCA -------------
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(object = filtered_seurat))
ElbowPlot(filtered_seurat)


# ----------- Clustering --------------
#Cluster cells that have similar expression pattern
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:15)
filtered_seurat <- FindClusters(filtered_seurat, resolution = c(0.1,0.3, 0.5, 0.8, 1))
View(filtered_seurat@meta.data)

DimPlot(filtered_seurat, group.by = "RNA_snn_res.0.5", label = TRUE)

# In this case---> set resolution to 0.5
# Set Idents to 0.5
Idents(filtered_seurat) <- "RNA_snn_res.0.5"



# -------- non-linear dimensional reduction ----------
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:15)
DimPlot(filtered_seurat, reduction = "umap")


# ----- Differential expression and marker selection  --------------
all.markers <- FindAllMarkers(filtered_seurat, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)

table(all.markers$cluster)

top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers

# -------- Cell annotation using singleR ----------
library(SingleR)
monaco.ref <- celldex::MonacoImmuneData()
sce <- as.SingleCellExperiment(DietSeurat(filtered_seurat))
sce

monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
table(monaco.main$pruned.labels)
table(monaco.fine$pruned.labels)

# Cell annotation seems to work becasue all are B cells (in the cellranger output only appeared B cells)

filtered_seurat@meta.data$monaco.main <- monaco.main$pruned.labels
filtered_seurat@meta.data$monaco.fine <- monaco.fine$pruned.labels

srat <- SetIdent(filtered_seurat, value = "monaco.fine")
DimPlot(srat, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()



##############################################################################################################
##############################################################################################################

############# BLOC 2. INTEGRATE GEX AND BCR #################

# Can be done with Immcantation, SCrepertoire ....

#--------------------------------------------------------------------

## Using IMMCANTATION: 
db_airr <- airr::read_rearrangement('/home/mvilardell/Documents/others/scellRNA/vdj/airr_rearrangement.tsv')

# In GEX: concatenate the sample ids with the cell ids (ensure uniqueness of cell barcodes across samples)
# In this particular case--> 1 sample 
#filtered_seurat[['cell_id_unique']]<-Cells(filtered_seurat) # or rownames of the metadata


# Examine GEX data:
filtered_seurat
names(filtered_seurat@commands)
names(filtered_seurat[[]]) # Be sure to have the cell_id_unique column


# examine a few columns of the BCR data
db_airr %>%
  ungroup %>%
  select(cell_id, v_call, j_call, c_call) %>%
  slice_sample(n = 5)

# make sure that there are no empty c_call rows
db_airr <- db_airr %>% dplyr::filter(!is.na(c_call), c_call != "")

# define custom colors for the cell types (for plotting later)
cols_anno <- c("Naive B cells" = "green", "Exhausted B cells" = "dodgerblue2", "Non-switched memory B cells" = "darkgoldenrod2",
               "Switched memory B cells" = "firebrick2")

metadata<-filtered_seurat@meta.data

#Idents(filtered_seurat)<- filtered_seurat$monaco.fine
filtered_seurat$annotated_clusters <- Idents(filtered_seurat)



# Check barcode styles
head(Cells(filtered_seurat))
head(db_airr$cell_id)

# Create new metadata columns -- match by using the unique cell ids
new_meta_cols <- data.frame(cell_id_unique = Cells(filtered_seurat))

# only add in the heavy chains (otherwise the row counts differ)
library(stringr)
bcr_data_selected <- db_airr %>%
  filter(str_detect(v_call, "^IGH"))

# select columns from BCR data
bcr_data_selected <- bcr_data_selected %>%
  select(cell_id, clone_id, c_call) 

colnames(bcr_data_selected)[1]<-'cell_id_unique'

# sort BCR data by cell id in GEX Seurat Object
new_meta_cols <- left_join(new_meta_cols, bcr_data_selected,
                           by = "cell_id_unique")

# integrate BCR data with the Seurat object
for (colname in colnames(bcr_data_selected)) {
  filtered_seurat[[colname]] <- new_meta_cols[[colname]]
}
filtered_seurat$contains_bcr <- !is.na(filtered_seurat$clone_id)


# examples of the new columns
ncol_meta <- ncol(filtered_seurat[[]])
filtered_seurat[[]][, (ncol_meta - 3):ncol_meta] %>%
  dplyr::filter(contains_bcr) %>% # let's look at the ones that match with BCRs
  slice_sample(n = 5) # row names = cell ids



# --- Examine integrated data ----
# How much of the GEX data has BCR
ggplot(filtered_seurat[[]], aes(x = contains_bcr, fill = c_call)) +
  geom_bar(color = "black", linewidth = 0.2) +
  labs(title = "GEX Data with Matching BCRs",
       x = "Contains BCRs", y = "Count")



# -------- Highligh BCR cells in the GEX UMAP --------
# have to set the idents
Idents(filtered_seurat) <- "annotated_clusters" # should be annotated

highlighted_cells <- Cells(filtered_seurat)[which(filtered_seurat$contains_bcr)]

UMAPPlot(object = filtered_seurat, cells.highlight = highlighted_cells,
         label = TRUE, label.size = 5, pt.size = 0.7) +
  scale_color_manual(name = "Data Type",
                     labels = c("non-BCR", "BCR"),
                     values = c("gray", "red"))


# ------- Integration of GEX cell annotations in the BCR data ----
# The annotation information of B cells (such as sub-types of B cells and 
# their associated UMAP coordinates) identified by the GEX data can be integrated into BCR data.

# select columns from the GEX data
gex_data_selected <-
  data.frame(cell_id_unique = filtered_seurat$cell_id_unique,
             gex_umap_1 = filtered_seurat@reductions$umap@cell.embeddings[, 1],
             gex_umap_2 = filtered_seurat@reductions$umap@cell.embeddings[, 2],
             gex_annotation = filtered_seurat$annotated_clusters)

# integrate GEX data with the BCR data
colnames(db_airr)[1]<-'cell_id_unique'
bcr_gex_data <- left_join(db_airr, gex_data_selected,
                          by = "cell_id_unique")

# keep BCR cells with matched GEX cells (assuming everything was annotated)
bcr_gex_data <- dplyr::filter(bcr_gex_data, !is.na(gex_annotation))



# Identify GEX clusters in the BCR UMAP
ggplot(data = bcr_gex_data %>% filter(str_detect(c_call, "^IGH")),
       aes(x = gex_umap_1, y = gex_umap_2, color = c_call)) +
  geom_point() +
  labs(x = "GEX UMAP 1", y = "GEX UMAP 2", color = "Isotype")


ggplot(data = bcr_gex_data %>% filter(str_detect(c_call, "^IGH")),
       aes(x = gex_umap_1, y = gex_umap_2, color = gex_annotation)) +
  geom_point() +
  labs(x = "GEX UMAP 1", y = "GEX UMAP 2", color = "Isotype")


# -----------------------------------------------------------------------------------------

# using SCRepertoire R package

filtered_contig_annotations <- read.csv("~/Documents/others/scellRNA/vdj/filtered_contig_annotations.csv")
combined.BCR <- combineBCR(BCR.contigs, 
                           samples = "P1", 
                           threshold = 0.85)



