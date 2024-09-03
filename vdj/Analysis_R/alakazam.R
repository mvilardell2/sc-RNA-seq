
################## TUTORIAL IMMCANTATION 

#######################################################################
#######################################################################

# Tutorial Alakazam

# Allows to investigate clonal lineages, diversity, gene usage and repertorire properties. 

library(airr)
library(alakazam)
data("ExampleDb")

# AIRR formatted files are loaded with the airr package:
db_airr <- airr::read_rearrangement('/home/mvilardell/Documents/others/scellRNA/vdj/airr_rearrangement.tsv')


# ------------- GENE USAGE ANALYSIS ---------------------
# Study the distribution of gene segments used in the recombination process of IG and TCRs. 
# Gene usage analysis examines which specific V, D, and J gene segments are selected during VDJ recombination.

# Important to have these columns: v_call, d_call, j_call

# Relative abundance and count of genes:
gene <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="gene")
gene

# Extract the family of the gene segment, order them and filter by a gene: 
ighv1 <- gene %>%
  mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
  filter(getFamily(gene) == "IGHV1")

# Plot V gene usage in the IGHV1 family by sample
g1 <- ggplot(ighv1, aes(x=gene, y=seq_freq)) +
  geom_point(aes(color=sample_id), size=5, alpha=0.8)


# Usage can also be quantified at the allele or family level (change the mode)
family <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="family")
family


# Quantify V family clonal usage by sample and isotype
family <- countGenes(ExampleDb, gene="v_call", groups=c("sample_id", "c_call"), 
                     clone="clone_id", mode="family")
head(family, n=4)
# Subset to IGHM and IGHG for plotting
family <- filter(family, c_call %in% c("IGHM", "IGHG"))
g3 <- ggplot(family, aes(x=gene, y=clone_freq)) +
  geom_point(aes(color=sample_id), size=5, alpha=0.8) +
  facet_grid(. ~ c_call)
plot(g3)



# ----------- Amino acid property analysis ---------------
# Analyse physicochemical properties of CDR3 regions is useful for several reasons: 
# - Infer how specific antibodies might bind to their respective antigens --> binding affinity and specificity
# - Depending on the structure (aa) of CDR3 loop, it can recognize a wide array of antigens. Antigens with similar sequences can be recognzed by the same sequence.
# - The properties of the region can help identifying targets, so can design peptides to modeulate the activity of these receptors
# - Also provide insights to the system's response, helps to design vaccines....

# Ammino acid properties
db <- ExampleDb[ExampleDb$sample_id == "+7d", ]# Get a subsample for faster analysis
db_props <- aminoAcidProperties(db, seq="junction", trim=TRUE, 
                                label="cdr3") #Adds new columns with the properties to the df
head(db_props)

# Generate plots: 
g1 <- ggplot(db_props, aes(x=c_call, y=cdr3_aa_length)) + 
  ggtitle("CDR3 length") + 
  xlab("Isotype") + ylab("Amino acids") +
  scale_fill_manual(name="Isotype", values=IG_COLORS) +
  geom_boxplot(aes(fill=c_call))
g1



# ---------------- Diversity analysis -------------------
# Variety and abundance of unique clones. Clone --> number of cells that correponds to the same clonotype (unique CDR3 region)
# A diversity repertoire-->respond to a wide range of pathogens.
# Changes in clonal diversity can indicate diseases.

# Clonal abundance: 
clones <- countClones(db_airr,group='c_call')
head(clones, 5)

# Clonal abundance distribution: 
curve <- estimateAbundance(ExampleDb, group="sample_id", ci=0.95, nboot=100, clone="clone_id")#calculates 95% interval
plot(curve, legend_title="Sample") # Rank abundance curve
# The most abundant clones are at the left of the x axis, and goes to the rigth with progessively less abundant clones. 
#  7d samples --> there is strong presence of the dominnt clones --> reduced diversity--> clonal response to a specific immune response
# -1h samples --> there are not many abundant clones -->more diversity

# Diversity curve: 
sample_curve <- alphaDiversity(ExampleDb, group="sample_id", clone="clone_id",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=100)
plot(sample_curve)#Diversity profile curve. Diversity across different diversity orders
# At q=0,the diversity index corresponds to species richness, where all clones are weighted equally regardless of their abundance.
# As q increases, the measure becomes more sensitive to the abundance of the most common clones, reducing the influence of rare ones.
# At q=0, red curve is higher --> higher species richness (more unique clones)
# As 1 increases, the curve dropes more in blue sample, so it has fewer dominant clones-->diversity driven by small number of abundant clones. 
#-1h samples--->higher diversity. 

# In this case, the -1h samples is more diverse. In 7d samples, some clones have become dominant --> diversity is reduced.



# -------------- Lineage reconstruction  --------------------
# Trace the evolutionary history by analyzing the sequences and genes of a clonotype.--->Process essential in B cells
# After B cells encounter an antigen -->somatic hypermutation happens in the V region --> increase the binding affinity.
# This process creates a family (clones), with slightly different antibody sequences. 
# B cells that producer high-affinity antibodies will survive and proliferate-->clonal expansion.

# Lineage tree-->evolutionary relationship between B cells in a clone


# Select a clone: 
sub_db <- subset(ExampleDb, clone_id == 3138)

# Preprocess clones (clean gaps, duplicated sequences ...)
clone <- makeChangeoClone(sub_db, text_fields=c("sample_id", "c_call"), 
                          num_fields="duplicate_count")
clone
# germline slot contains the inferred ancestral sequence before SHM

clone@data[,c('sample_id','c_call','duplicate_count')]

# Run PHYLIP and parse output
# Need PHYLIP APP
phylip_exec <- "~/apps/phylip-3.69/dnapars"
graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)
# ....

