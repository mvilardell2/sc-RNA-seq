###########################################################################################################
############################################################################################################

# Tutorial scRepertoire

# Single cell immune profiling


library(scRepertoire)

# Input--> filtered_contig_annotations.csv

###### LOAD DATA ########
S1 <- read.csv("filtered_contig_annotations.csv")
head(S1)
contig_list2 <- list(S1)
names(contig_list2)<-'S1'

# Example data: 
data('contig_list')
head(contig_list[[1]])


####### Combine contigs into clones  #########

table(contig_list2[[1]]$chain)

# - TCR data --> TRA,TRB,TRG ...
# - BCR data --> IGH, IGK, IGL
# My data only contain BCR data

combined.BCR <- combineBCR(contig_list2, 
                           samples = 'S1',
                           threshold = 0.85)

head(combined.BCR[[1]])

#remove S1_ prefix 
combined.BCR$S1$barcode<-gsub('S1_','',combined.BCR$S1$barcode)



########## BASIC CLONAL VISUALIZATIONS #################

# Explore the clones. Return the relative percent or total unique clones.
clonalQuant(combined.BCR, 
            cloneCall="strict", 
            chain = "both", # IGH or IGL to select a single chain 
            scale = FALSE)



########### COMBINE GEX + BCR DATA ############
# Cell barcodes of the metadata of the GEX and BCR data must match. 
sce <- combineExpression(combined.BCR, 
                         filtered_seurat, 
                         cloneCall="gene")

head(sce)

DimPlot(sce, reduction = "umap", group.by = "c_call")

b_cell_markers <- c("CD79A","CD79B")
FeaturePlot(sce, features = b_cell_markers)


# Find markers associated witht the Clonotypes that contain IGHV1
sce$cells_of_interest <- FALSE
sce$cells_of_interest[grep("IGHV1", sce$CTstrict)] <- TRUE

head(sce@meta.data)
table(sce$cells_of_interest)
Idents(sce) <- sce$cells_of_interest
DimPlot(sce)

FM <-FindMarkers(sce, ident.1 = "TRUE")


seurat <- highlightClones(sce, cloneCall= "aa", 
                          sequence = c("CARGDSSGWRGGNWFDPW_CQSYDSSLSDVF", "CAMGYCINNNCYEGWFDPW_CQQYYDTPRTF"))
DimPlot(seurat, group.by = "highlight")
