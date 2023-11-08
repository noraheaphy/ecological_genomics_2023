# load packages
library(WGCNA)
options(stringsAsFactors = FALSE)
library(DESeq2)
library(ggplot2)
library(tidyverse)

# install.packages("remotes")
# remotes::install_github("kevinblighe/CorLevelPlot")

library(CorLevelPlot) 
library(gridExtra)
library(Rmisc) 
library(pheatmap)

############################## Day 1 - 10/25/23 ##############################

##### 1. Import the counts matrix and metadata and filter using DESeq2 #####

countsTable <- read.table("myresults/salmon.isoform.counts.matrix.filteredAssembly", 
                          header = TRUE, row.names = 1)
head(countsTable)
dim(countsTable)

# bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
countsTableRound <- round(countsTable) 
head(countsTableRound)

#import the sample description table
# conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
# head(conds)

sample_metadata = read.table(file = "mydata/Ahud_trait_data.txt",header = T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = sample_metadata, 
                              design = ~ 1)

dim(dds)

# Filter out genes with too few reads:
# remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

##### 2. QC - outlier detection #####

# detect outlier genes
gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples) # any that are false here, we would need to filter out

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree) 

# pca - method 2
pca <- prcomp(t(countsTable))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

##### 3. Normalization #####

colData <- row.names(sample_metadata)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(countsTableRound)) # to see if all samples are present in both
all(rownames(colData) == colnames(countsTableRound))  # to see if all samples are in the same order



# perform variance stabilization
dds_norm <- vst(dds)
# dds_norm <- vst(normalized_counts)

# get normalized counts
norm.counts <- assay(dds_norm) %>% t()


##### 4. Network Construction #####

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function; this step takes a couple minutes
# signed means cares about whether downregulated vs. upregulated
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

# Warning message:
# executing %dopar% sequentially: no parallel backend registered 

sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)

# based on this plot, choose a soft power to maximize R^2 (above 0.8) and minimize connectivity
# for these ahud data: 6-8; Higher R2 should yield more modules.


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 6
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package


# this step also takes a few minutes; ideally your maxBlockSize is larger than your number of 
# genes to run the memory-intensive network construction all at once.
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 26000, # set lower to save computational power, computer will run in chunks
                          minModuleSize = 30, # to call something a module, needs >=30 genes
                          reassignThreshold = 0,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)

# TOMtype (Topological Overlap Matrix type) parameter - unsigned - doesn't consider positive/negative co-expression
# signed - when you want to consider the direction of co-expression interaction, e.g., activating or inhibiting
# WGCNA often uses a dendrogram-based approach to identify modules. The choice of the 
# height cut in the dendrogram can determine the number of modules. Selecting a higher
# cut height results in fewer, larger modules, while a lower cut height leads to more, 
# smaller modules.

cor <- temp_cor

############################## Day 2 - 10/30/23 ##############################

# To load the object from last time
bwnet <- readRDS("myresults/bwnet.rds")

##### 5. Model eigengenes #####

module_eigengenes <- bwnet$MEs
head(module_eigengenes)

# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

# grey module = all genes that don't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module.

##### 6. Relate modules to traits #####

# module trait associations
traits <- sample_metadata[, c(5,8,11,14,17)]


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# test for correlation (Pearson's) between eigengenes and trait data
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

names(heatmap.data)

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[12:16],
             y = names(heatmap.data)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"))


module.gene.mapping <- as.data.frame(bwnet$colors) # assigns module membership to each gene
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames()

groups <- sample_metadata[,c(3,1)]
module_eigengene.metadata <- merge(groups, heatmap.data, by = 'row.names')

# create a summary data frame of a particular module eigengene information
MEyellow_summary <- summarySE(module_eigengene.metadata, measurevar="MEyellow", 
                              groupvars = c("Generation","treatment"))

#Plot a line interaction plot of a particular module eigengene
ggplot(MEyellow_summary, aes(x=as.factor(Generation), y=MEyellow, color=treatment, fill = treatment, shape = treatment)) +
  geom_point(size=5, stroke = 1.5 ) +
  geom_errorbar(aes(ymin=MEyellow-se, ymax=MEyellow+se), width=.15) +
  geom_line(aes(color=treatment, group=treatment, linetype = treatment)) +
  scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  xlab("Generation") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4)) +
  theme(text = element_text(size = 20)) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6))

# 6B. Intramodular analysis: Identifying driver genes ---------------

# Get top hub genes (genes with highest connectivity in the network)
hubs <- chooseTopHubInEachModule(norm.counts, bwnet$colors, type = "signed", omitColors = "")
hubs

### Plot Individual genes  to check! ### 

d <-plotCounts(dds, gene = "TRINITY_DN11845_c0_g1::TRINITY_DN11845_c0_g1_i9::g.36434::m.36434", intgroup = (c("treatment","Generation")), returnData = TRUE)
d_summary <- summarySE(d, measurevar = "count", groupvars=c("Generation","treatment"))

# here y-axis is normalized number of reads that map to the particular gene we're looking at
# that's correlated strongly with this particular module.

ggplot(d_summary, aes(x = as.factor(Generation), y=count, color=treatment, fill = treatment, shape = treatment)) +
  geom_point(size=5, stroke = 1.5 ) +
  geom_errorbar(aes(ymin=count-se, ymax = count+se), width=.15) +
  geom_line(aes(color=treatment, group = treatment, linetype = treatment)) +
  scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) +
  scale_shape_manual(values = c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  xlab("Generation") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6))

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10,1:10]

## Now with purple

# Make a heat map of gene expressions within modules.
# Use the norm.counts matrix, subset based on module membership
t_norm.counts <- norm.counts %>% t() %>% as.data.frame()

# Purple module
purple_transcripts <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'purple') %>% 
  rownames()

t_norm.counts_purple <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% purple_transcripts)

t_norm.counts_purple <- t_norm.counts_purple - rowMeans(t_norm.counts_purple)
df <- as.data.frame(colData(dds)[,c("Generation","treatment")])

#blue to yellow color scheme
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "yellow"))(paletteLength)
myBreaks <- c(seq(min(t_norm.counts_purple), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(t_norm.counts_purple)/paletteLength, max(t_norm.counts_purple), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_purple, color = myColor, breaks = myBreaks,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Yellow")
