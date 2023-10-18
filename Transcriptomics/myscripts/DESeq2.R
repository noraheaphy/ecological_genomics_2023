if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("hexbin")
BiocManager::install("vsn")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)

###################### Import Data ######################

# Import the counts matrix
countsTable <- read.table("myresults/salmon.isoform.counts.matrix", header = TRUE, row.names = 1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample discription table
conds <- read.delim("myresults/ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

# Let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))

###################### Explore Data Counts ######################

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound),
        cex.names = 0.5, las = 3, ylim = c(0,21000000))
abline(h = mean(colSums(countsTableRound)), col = "blue", lwd = 2)

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # [1] 11930.81 - tonsa, 6076.078 - hudsonica genes, 2269 - hudsonica isoform
median(rowSums(countsTableRound)) # [1] 2226 - tonsa, 582 - hudsonica, 109

apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound,1,mean),xlim=c(0,1000), ylim=c(0,120000),breaks=10000)

#### Create a DESeq object and define the experimental design here with the tilda ####

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ generation + treatment)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 30) >= 28,]
nrow(dds) 

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

