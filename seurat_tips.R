#########  SEURAT TIPS and COMMANDS  ###########
################################################

library(Seurat)
library(SeuratData)
InstallData("pbmc3k")
pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")

set.seed(42)
pbmc$replicate <- sample(c("rep1", "rep2"), size = ncol(pbmc), replace = TRUE)


################# SEURAT OBJECT ##################

# SLOTS: 
# - assays
# - meta.data	
# - active.assay:	Name of active, or default, assay
# - active.ident:	Identity classes for the current object
# - graphs:	A list of nearest neighbor graphs
# - reductions:	A list of DimReduc objects
# - project.name:	User-defined project name (optional)
# - tools	Empty list. Tool developers can store any internal data from their methods here
# - misc	Empty slot. User can store additional information here



# --- Basic manipulation ---
str(pbmc)

# Get cell names
colnames(pbmc)
Cells(pbmc)

# Get features names: 
Features(pbmc)
rownames(pbmc)

# Number of features: 
nrow(pbmc)

# Number of cells: 
ncol(pbmc)
dim(pbmc)

# Head of rows and cells
head(rownames(pbmc))
head(colnames(pbmc))

# Names associated with objects
names(pbmc)
# This can be passed to the [[ ]] to extract them.


# list of object layers (within the RNA assay)
Layers(pbmc)

# working with multimodal objects list assays
Assays(pbmc)



# ---- Data access ----
# Pulling specific Assay, DimReduc, or Graph objects can be done with [[ ]]
pbmc[['RNA']]
pbmc[['pca']]

# Retrieve data in an expression matrix RNA counts matrix
pbmc[["RNA"]]$counts

# Access data from counts, data, or scale.data slots
GetAssayData(object = pbmc, slot = 'scale.data')[1:3, 1:3]

# Get metadata----> can acces with @ and $ or with [[ ]]
pbmc@meta.data
pbmc[[ ]]
head(pbmc[['nCount_RNA']])
head(pbmc$nCount_RNA)

pbmc@meta.data$orig.ident


# HVFInfo pulls feature mean and disperson from an ASSAY
#Useful for viewing the results of FindVariableFeatures
head(x = HVFInfo(object = pbmc))



###############################################################
##################  ASSAY DATA  ###############################


# Stores single cell data
# It has a single assay (in a basic experiment) --> RNA
# This assay stores multiple transformation of the data (counts, data, and scale.data)


# SLOTS: 
# - counts
# - data
# - scale.data
# - Key
# - var.features

rna<- pbmc[['RNA']]
dim(rna)
str(rna)

# ------ Expression data access ------
# Pull expression data from the data slot: 
rna@layers$data[1:3,1:3]
# or
pbmc[["RNA"]]$counts
# or 
rna$counts

# Alternate accessor function with the same result
LayerData(pbmc, assay = "RNA", layer = "counts")

# Expression data is accesed with GetAssayData()
GetAssayData(object = rna, slot = 'scale.data')[1:3, 1:3]

## Set assay data: 
# Set expression data assume new.data is a new expression matrix
pbmc[["RNA"]]$counts <- new.data

# SetAssayData from Seurat v4 is still supported
pbmc <- SetAssayData(object = pbmc, slot = "counts", new.data = new.data)



# --------- Subset Seurat Object -----------

# Subset Seurat object based on identity class, also see ?SubsetData
subset(x = pbmc, idents = "B")

# Subset on the expression level of a gene/feature
subset(x = pbmc, subset = MS4A1 > 2.5)

# Subset on a combination of criteria
subset(x = pbmc, subset = MS4A1 > 2.5 & PC_1 > 5)


# --------------------------------------

# Feature-level meta data ---> [[ ]]
colnames(rna[[ ]])

# Get the vector of variable features
VariableFeatures(pbmc)
head(x = VariableFeatures(object = rna))

# KEY: 
# Access and sets the key slot for the assay
Key(rna)
Key(object = rna) <- 'myRNA_'
Key(object = rna)




# ----- Identity class labels ----------

# Setting and retrieving cell identities

# Set identity classes to an existing column in meta data
Idents(object = pbmc) <- "seurat_annotations"

# View cell identities, get summary table
Idents(pbmc)
table(Idents(pbmc))

# Set identity to CD4 T cells for all cells
Idents(pbmc) <- "CD4 T cells"

# Set for a selected group of cells
pbmc.cells <- Cells(pbmc)
Idents(object = pbmc, cells = pbmc.cells[1:10]) <- "CD4 T cells"

# Get cell identity classes
Idents(object = pbmc)
levels(x = pbmc)

# Put cell identity classes in metadata
pbmc[["old.ident"]] <- Idents(object = pbmc)

# Rename identity classes
pbmc <- RenameIdents(object = pbmc, `CD4 T cells` = "T Helper cells")



# -------- See the methods available ----------
utils::methods(class = 'Assay')



#################################################
############  DimReduc slot  ###################

# This slot represents a dimensional reduction object

str(pbmc)
pca<-pbmc[['pca']] # This is found within pbmc@reductions

length(pca)# Number of dimensions

# SLOTS: 
# - cell.embeddings:	A matrix with cell embeddings
# - feature.loadings:	A matrix with feature loadings
# - feature.loadings.projected:	A matrix with projected feature loadings
# - assay.used:	Assay used to calculate this dimensional reduction
# - stdev:	Standard deviation for the dimensional reduction
# - key:	A character string to facilitate looking up features from a specific DimReduc
# - jackstraw:	Results from the JackStraw function


# ----- Data access -----

# Pulling features can be done: 
pca[1:3,1:3]

# Features can also be accessed with Loadings()
Loadings(object = pca, projected = FALSE)[1:3, 1:3]

# Cell Embeddings: 
pca[[1:3, 1:3]]
Embeddings(object = pca)[1:3, 1:3]

# Key
Key(object = pca)

# Stdev gets the vector of standard deviations for each dimension embedded.
Stdev(object = pca)



#########################################################
#########################################################

# -------- Switch identity class -------------
# Plot UMAP, coloring cells by cell type (currently stored in object@ident)
DimPlot(pbmc, reduction = "umap")

##  UMAP colored by replicated
pbmc$CellType <- Idents(pbmc)
head(pbmc@meta.data)
# switch the identity class of all cells to reflect replicate ID
Idents(pbmc) <- "replicate"
DimPlot(pbmc, reduction = "umap")


# Switch back to cell type labels
Idents(pbmc) <- "CellType"

# Alternatively: 
DimPlot(pbmc, reduction = 'umap', group.by = 'replicate') 


# ------- Tabulate cells by cluster ID, replicate, or both ------

# How many cells are in each replicate?
table(pbmc$replicate)

# What proportion of cells are in each cluster?
prop.table(table(Idents(pbmc)))

# How does cluster membership vary by replicate?
table(Idents(pbmc), pbmc$replicate)


# ---------- Selecting cells and subsetting -----------

# What are the cell names of all NK cells?
WhichCells(pbmc, idents = "NK")


# How can I extract expression matrix for all NK cells (perhaps, to load into another package)
nk.raw.data <- as.matrix(GetAssayData(pbmc, slot = "counts")[, WhichCells(pbmc, ident = "NK")])

# Can I create a Seurat object based on expression of a feature or value in object metadata?
subset(pbmc, subset = MS4A1 > 1)


# Can I create a Seurat object of just the NK cells and B cells?
subset(pbmc, idents = c("NK", "B"))


# -------- Calculating the average gene expression within a cluster-----------

cluster.averages <- AverageExpression(pbmc)
head(cluster.averages[["RNA"]][, 1:5])

## Return this information as a Seurat object
# replace spaces with underscores '_' so ggplot2 doesn't fail
orig.levels <- levels(pbmc)
Idents(pbmc) <- gsub(pattern = " ", replacement = "_", x = Idents(pbmc))


cluster.averages <- AverageExpression(pbmc, return.seurat = TRUE)
cluster.averages

DoHeatmap(cluster.averages, features = unlist(TopFeatures(pbmc[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)
