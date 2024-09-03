#################################################################################################################
#################################################################################################################
# Tutorial immunarch

library(immunarch)
immdata_10x <- repLoad('/home/mvilardell/Documents/others/scellRNA/vdj/')

names(immdata_10x$data)<-'sample1'
df<-data.frame(Sample='sample1',Sex='M')
immdata_10x$meta<-df

# I will use the tutorial data for better understanding

data('immdata')

############ Exploratory analysis ##############

# Get the chains of the receptor of each cell
first_part <- unlist(strsplit(immdata_10x$data$sample1$J.name, ";"))
result <- substr(sub(";.*", "", first_part), 1, 3)
table(result)


# Number of clonotypes + statistical tests

exp_vol <- repExplore(immdata$data, .method = "volume") #total number of unique clonotypes
p1 <- vis(exp_vol, .by = c("Status"), .meta = immdata$meta)
p2 <- vis(exp_vol, .by = c("Status", "Sex"), .meta = immdata$meta)
p1 + p2

# Comparing number of reads / UMIs 

exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_cnt <- repExplore(immdata$data, .method = "count")
exp_vol <- repExplore(immdata$data, .method = "volume")

p1 <- vis(exp_len)
p2 <- vis(exp_cnt)
p3 <- vis(exp_vol)

p1

p2+p3


# ------------- Clonality ------------

# Estimate diversity of samples. How much the repertoire is dominate by few clones. Measure the amount of the most or least frequent clonotypes.
imm_pr <- repClonality(immdata_10x$data, .method = "clonal.prop")
imm_pr
# How much the repertoire is occupied by the most frequent clones
# MS1-->only 2 clonotypes in the top 10% of the clones. Suggest that has a few highly frequent clones.
# Samples with low number of clones (or low proportion) suggest that it is highly clonal, so few clones dominate the repertoire. 
# Samples with high number of clones suggest that it has a more diverse response, due to 10% of the abundance are made up of a large amount of clonotypes.


# Calculate the proportion of the repertoire occupied by the most abundant clones at diff tresholds-->study the dominance of the most abundant
imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top
# Column 10 shows the proportion of the repertoire occupied by the top 10 most abundant clones
# High proportion (MS1) -->indicates a more clonal response, few clones dominate. So the 10 most abundant clones make the 20% of the repertoire (MS1)
# Small number of clones dominate the repertoire (MS1).

vis(imm_top) + vis(imm_top, .by = "Status", .meta = immdata$meta)




############ Repertoire overlap and public clonotypes ######################
# To measure repertoire similarity between samples

imm_ov1 <- repOverlap(immdata$data, .method = "public", .verbose = F)
imm_ov2 <- repOverlap(immdata$data, .method = "morisita", .verbose = F)

p1 <- vis(imm_ov1)
p1
# Each entry represents the number of shared clonotypes between the samples. Higher number indicate similar repertoires. 

vis(imm_ov1, "heatmap2") # can also be ploted as a heatmap

p2 <- vis(imm_ov2, .text.size = 2)
p2
# It returns a value that indicates the similarity between the samples. This is based on their clonotypic distribution. 
# Higher values, higher similarity. 


## Analyse the overlap measures: 
# Apply different analysis algorithms to the matrix of public clonotypes:
# "mds" - Multi-dimensional Scaling--reduces the high-dimensional overlap to 2 dimensions.
repOverlapAnalysis(imm_ov1, "mds")

# Visualise the results
repOverlapAnalysis(imm_ov1, "mds") %>% vis()

# Samples with similar overlap will be closer toghether, while those with different profiles will be further apart. 


# Clusterise the MDS resulting components using K-means
repOverlapAnalysis(imm_ov1, "mds+kmeans") %>% vis()



# Build a table will all clonotypes
pr.nt <- pubRep(immdata$data, "nt", .verbose = F) # nt --> build the table using CDR3 nucleotide sequences




########### Gene Usage Analysis ###########
# Determine with gene segment is present in each VDJ-recombination-segment of each cell

gene_stats()#Data table containing gene segments that are associated with IG and TCR for several species.
# For example, in Homo sapiens, there are 30 distincs D gene segments for the heavy chain, 13 J gene segments and 248 V gene segments...
# H---> heavy chain
# K; kappa, L; lambda ligth chain 

imm_gu <- geneUsage(immdata_10x$data, "hs.trbv")
# Retrieves the amount of TRBV genes in the dataset. Higher counts-->the gene is more frequently used.
# This can reveal which genes are more commonly used in each sample. This is important for understanding cell diversity. 

# Compute the distribution of the first two samples
imm_gu$Names<-sapply(strsplit(strings, ";"), `[`, 1)
imm_gu <- geneUsage(immdata_10x$data[c(1, 2)], "hs.trbv", .norm = T)
vis(imm_gu)

imm_gu <- geneUsage(immdata$data, "hs.trbv", .norm = T)
vis(imm_gu, .grid = T)




############ Diversity estimation ###########
# Estimation of the repertoire diversity

# Compute statistics and visualise them
# Chao1 diversity measure
div_chao <- repDiversity(immdata$data, "chao1") #Chao1 is an estimator for species richness (unique clonotypes)
# - Estimator: number of unique clonotypes (estimation)

p1 <- vis(div_chao)
p2 <- vis(div_chao, .by = c("Status", "Sex"), .meta = immdata$meta)

p1
p2

# Hill numbers
div_hill <- repDiversity(immdata$data, "hill") #Quantify diversity considering oth tichness and evennes (distribution of clonotype freq)

p3 <- vis(div_hill, .by = c("Status", "Sex"), .meta = immdata$meta)
p3


# D50
div_d50 <- repDiversity(immdata$data, "d50")
# This indicate how many unique clones contribute to 50 % of the diversity. Higher values-->more diversity, meaning it take more unique clones to account for the 50%
# Lower values--> less diversity, fewer clones dominate the repertoire.
p4 <- vis(div_d50)
p5 <- vis(div_d50, .by = "Status", .meta = immdata$meta)


# Ecological diversity measure
div_div <- repDiversity(immdata$data, "div")
p6 <- vis(div_div)
p6
# Higher values grater diversity, immune repertoire is more evenly distributed among many different clones. 
# Lower values, less diversity, indicates a more focused immune response (example: particular infection or immune condition)


# Rarefaction analysis --> to understand how the diversity would change if different number of sequences were present in the sample. 
imm_raref <- repDiversity(immdata$data, "raref", .verbose = F)

p1 <- vis(imm_raref)
p2 <- vis(imm_raref, .by = "Status", .meta = immdata$meta)

p1 + p2


repDiversity(immdata$data, "raref", .verbose = F) %>% vis(.log = TRUE)




############## Tracking clonotypes ######################
# Track the changes in the frequency of clonotypes. Used to understand how specific clones respond to infections...
# Example: study how a clonotype behaves in response to a pathogen or in a disease. 

# ------ Tracking the most abundance clonotypes -------
# Choose the top 5 most abundant clonotypes from the first repertoire (sample), and track them using the CDR3 nucleotide sequence: 
tc1 <- trackClonotypes(immdata$data, list(1, 5), .col = "nt")
p1 <- vis(tc1)

p1

# ----- Tracking clonotypes with specific nucleotide or aa -----
target <- c("CASSLEETQYF", "CASSDSSGGANEQFF", "CASSDSSGSTDTQYF", "CASSLAGGYNEQFF", "CASSDSAGGTDTQYF", "CASSLDSYEQYF", "CASSSAGGYNEQFF")
tc <- trackClonotypes(immdata$data, target, .col = "aa")
vis(tc)

# ---- Tracking clonotypes with specific sequences and genes ----
# Choose the 10 most abundant clonotypes and its gene segments from the first repertoire
target <- immdata$data[[1]] %>%
  select(CDR3.aa, V.name) %>%
  head(10)

target

tc <- trackClonotypes(immdata$data, target)
vis(tc)


# Add information about timepoints: 
immdata$meta$Timepoint <- sample(1:length(immdata$data))
immdata$meta
sample_order <- order(immdata$meta$Timepoint)
target <- c("CASSLEETQYF", "CASSDSSGGANEQFF", "CASSDSSGSTDTQYF", "CASSLAGGYNEQFF", "CASSDSAGGTDTQYF", "CASSLDSYEQYF", "CASSSAGGYNEQFF")
tc <- trackClonotypes(immdata$data, target, .col = "aa")
vis(tc, .order = sample_order)




############# Annotate clonotypes ####################
# AIRR databases: VDJDB, McPAS-TCR, PIRD TBAdb

# VDJBD: database for TCR sequences
vdjdb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb", .species = "HomoSapiens", .chain = "TRB", .pathology = "CMV")
vdjdb

mcpas = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/McPAS-TCR.csv.gz", "mcpas", .species = "Human", .chain = "TRB", .pathology = "Cytomegalovirus (CMV)")
mcpas

dbAnnotate(immdata$data, vdjdb, "CDR3.aa", "cdr3")

dbAnnotate(immdata$data, mcpas, c("CDR3.aa", "V.name"), c("CDR3.beta.aa", "TRBV"))

# Create a database with specific set of sequences: 
local_db = data.frame(Seq = c("CASSDSSGGANEQFF", "CSARLAGGQETQYF"), V = c("TRBV6-4", "TRBV20-1"), stringsAsFactors = F)
dbAnnotate(immdata$data, local_db, c("CDR3.aa", "V.name"), c("Seq", "V"))




############## Kmer and sequence motif analysis ##################
# Kmers--> sequences of leng K that are extracted from a longer sequence

# ----------- Kmer -------------
# For one sample: 
kmers <- getKmers(immdata$data[[1]], 3)
kmers

# For a batch of repertoires: 
kmers <- getKmers(immdata$data, 5)
kmers
vis(kmers)

p1 <- vis(kmers, .head = 5)
p1
p3 <- vis(kmers, .head = 10, .position = "dodge")
p3


# ------------ Sequence motifs analysis ---------
kmers <- getKmers(immdata$data[[1]], 5)
kmer_profile(kmers)
# Indicate how frequently each kmer occuers in different positions.
# Ex: A appears 8955 time in first position.
# This can be used to discover motifs that are overrepresented
kp <- kmer_profile(kmers, "self")
p1 <- vis(kp)
p2 <- vis(kp, .plot = "seq")

p1 + p2







####################### Clustering ##########################
# 1. Calculate distances between sequences
# 2. Clustering

data("bcrdata")

# ----------- Calculating distances ----------------
# CDR3 is the most variable region in BCR and TCR, and it is the region used to calculate distances. 

TCRdata <- repFilter(immdata, .method = "by.meta", .query = list(Sample = include("A2-i129")))$data
BCRdata <- bcrdata$data

#using aa CDR3 sequences for calculating distance
distTCR <- seqDist( TCRdata, .col = 'CDR3.aa')

#using nt CDR3 sequences for calculating distance
distTCR <- seqDist( TCRdata, .col = 'CDR3.nt')

#using full sequences for calculating distance
distBCR <- seqDist( BCRdata, .col = 'Sequence')


# ----------------- Clustering ------------------
# Find TCR and BCR that have similar CDR3 sequences. 

#clustering TCR by CDR3 regions
clustTCR <- seqCluster(TCRdata, distTCR, .perc_similarity = 0.9)

#calculate number of unique clusters in column 'Cluster'
clustTCR$"A2-i129" %>% .$Cluster %>% unique() %>% length()






###################### BCR pipeline #######################
data("bcrdata")

# --------- Reconstructing clonal lineages -------------
# Clonal lineages represent B cells that have a common origin, so they are clustered according to their similarity.
# 1. Calculate distances. 2. clustering. 

#calulate distance matrix
distBCR <- seqDist(bcrdata$data %>% top(500))

#find clusters
bcrdata$data <- seqCluster(bcrdata$data %>% top(500), distBCR, .perc_similarity = 0.9)



# --------- Building germline --------------
# Each lineage has its own germline sequence that represents the ancestral sequence (sequence after VDJ recombination but before
# B cell maturarion process)
# It is useful for assesing the degree of mutation and maturity: 


#germline example
bcrdata$data %>%
  top(1) %>%
  repGermline(.threads = 1) %>% .$full_clones %>% .$Germline.sequence
# A germline is represented by V gen -N sequence and J gene

