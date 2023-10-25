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
