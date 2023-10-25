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
library(pheatmap)
library(RColorBrewer)
library(eulerr)

###################### Import data ######################

# Import the counts matrix
countsTable <- read.table("myresults/salmon.isoform.counts.matrix.filteredAssembly", 
                          header = TRUE, row.names = 1)
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

###################### Explore data counts ######################

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound),
        cex.names = 0.5, las = 3, ylim = c(0,21000000))
abline(h = mean(colSums(countsTableRound)), col = "blue", lwd = 2)

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # 8217.81
median(rowSums(countsTableRound)) # 377

apply(countsTableRound, 2, mean) # 2 in the apply function does the action across columns
apply(countsTableRound, 1, mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound, 1, mean), xlim = c(0,1000), ylim = c(0,60000), breaks = 1000)

# Create a DESeq object and define the experimental design here with the tilda

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ generation + treatment)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

###################### Visualize data quality ######################

# Check the quality of the data by sample clustering and visualization
# The goal of transformation "is to remove the dependence of the variance on the mean, 
# particularly the high variance of the logarithm of count data when the mean is low."

# this gives log2(n + 1)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
meanSdPlot(assay(vsd))

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$generation, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# AM_F11_Rep3 and OA_F2_Rep2 could be outliers
# OA = ocean acidification
# OW = ocean warming
# OWA = ocean warming and acidification
# AM = ambient control

###################### Visualize global gene expression with PCA ######################

# PCA to visualize global gene expression patterns

# first transform the data for plotting using variance stabilization
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = c("treatment", "generation"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = treatment, shape = generation)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar[2],"% variance")) + 
  coord_fixed()

# Let's plot the PCA by generation in four panels

data <- plotPCA(vsd, intgroup = c("treatment", "generation"), returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

dataF0 <- subset(data, generation == 'F0')

F0 <- ggplot(dataF0, aes(PC1, PC2)) +
  geom_point(size = 10, stroke = 1.5, aes(fill = treatment, shape = treatment)) +
  xlab(paste0("PC1: ", percentVar[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar[2],"% variance")) +
  ylim(-10, 25) + xlim(-40, 10)+ # zoom for F0 with new assembly
  #ylim(-40, 25) + xlim(-50, 50)+ # new assembly limits
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values = c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", 
                                                         "OWA"))+
  scale_fill_manual(values = c('#6699CC', "#F2AD00", "#00A08A", "#CC3333"), 
                    labels = c("Ambient", "Acidification","Warming", "OWA"))+
  ##theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

F0


#png("PCA_F0.png", res=300, height=5, width=5, units="in")
#ggarrange(F0, nrow = 1, ncol=1)
#dev.off()

#### F2

dataF2 <- subset(data, generation == 'F2')

F2 <- ggplot(dataF2, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill = treatment, shape = treatment)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23), labels = c("Ambient", "Acidification","Warming"))+
  # scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) + 
  #scale_color_manual(values=c('black')) +
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A"), 
                    labels = c("Ambient", "Acidification","Warming"))+
  theme(legend.position = c(0.83, 0.85), legend.background = element_blank(), 
        legend.box.background = element_rect(colour = "black")) +
  #scale_size(guide="none") +
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

F2


# png("PCA_F2.png", res=300, height=5, width=5, units="in")
# ggarrange(F2, nrow = 1, ncol=1)
# dev.off()

# Yes - F2 is missing one ambient replicate

#### F4

dataF4 <- subset(data, generation == 'F4')

F4 <- ggplot(dataF4, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+ # limits with filtered assembly
  #ylim(-20, 10) + xlim(-40, 25)+  # zoom with filtered assembly
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  # scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) + 
  #scale_color_manual(values=c('black')) +
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  #theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #scale_size(guide="none") +
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

F4


# png("PCA_F4.png", res=300, height=5, width=5, units="in")
# ggarrange(F4, nrow = 1, ncol=1)
# dev.off()


#### F11

dataF11 <- subset(data, generation == 'F11')

F11 <- ggplot(dataF11, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-45, 25) + xlim(-50, 55)+
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,24), labels = c("Ambient", "OWA"))+
  scale_fill_manual(values=c('#6699CC', "#CC3333"), labels = c("Ambient", "OWA"))+
  guides(shape = guide_legend(override.aes = list(shape = c( 21, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

F11


# png("PCA_F11.png", res=300, height=5, width=5, units="in")
# ggarrange(F11, nrow = 1, ncol=1)
# dev.off()

ggarrange(F0, F2, F4, F11, nrow = 2, ncol = 2)


###################### Explore differential expression ######################

## Check on the DE results from the DESeq command way above ##

resAM_OWA <- results(dds, name = "treatment_OWA_vs_AM", alpha = 0.05)

resAM_OWA <- resAM_OWA[order(resAM_OWA$padj),]
head(resAM_OWA)  

summary(resAM_OWA)


resAM_OW <- results(dds, name = "treatment_OW_vs_AM", alpha = 0.05)

resAM_OW <- resAM_OW[order(resAM_OW$padj),]
head(resAM_OW)  

summary(resAM_OW)

### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene = "TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.747::m.747", 
               intgroup = (c("treatment", "generation")), returnData = TRUE)
d

p <-ggplot(d, aes(x = treatment, y = count, color = treatment, shape = generation)) + 
  theme_minimal() + theme(text = element_text(size=20), 
                          panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position = position_jitter(w = 0.2, h = 0), size = 3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size = 5, alpha = 0.7) 
p

###################### Run another model ######################

#################### MODEL NUMBER 2 - subset to focus on effect of treatment for each generation

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design= ~ treatment)

dim(dds)

# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_sub <- subset(dds, select = generation == 'F0')
dim(dds_sub)

# Perform DESeq2 analysis on the subset
dds_sub <- DESeq(dds_sub)

resultsNames(dds_sub)

res_F0_OWvAM <- results(dds_sub, name = "treatment_OW_vs_AM", alpha = 0.05)

res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM)

summary(res_F0_OWvAM)


### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <- plotCounts(dds_sub, gene = "TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130", 
               intgroup = (c("treatment","generation")), returnData = TRUE)
d

p <- ggplot(d, aes(x = treatment, y = count, color = treatment, shape = generation)) + 
  theme_minimal() + theme(text = element_text(size=20), 
                          panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position = position_jitter(w = 0.2, h = 0), size = 3)
p <- p + stat_summary(fun = mean, geom = "point", size = 5, alpha = 0.7) 
p

# We can make an MA plot
plotMA(res_F0_OWvAM, ylim = c(-4,4))

# Heatmap of top 20 DEGs sorted by pvalue

# By environment
vsd <- vst(dds_sub, blind = FALSE)

topgenes <- head(rownames(res_F0_OWvAM),100)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_sub)[,c("generation","treatment")])
pheatmap(mat, annotation_col = df)
pheatmap(mat, annotation_col = df, cluster_cols = F)


###################### Plot overlapping DEGs in a Venn Euler diagram ######################

# For OW vs AM
res_F0_OWvAM <- results(dds_sub, name="treatment_OW_vs_AM", alpha=0.05)
res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM)

summary(res_F0_OWvAM)
res_F0_OWvAM <- res_F0_OWvAM[!is.na(res_F0_OWvAM$padj),]
degs_F0_OWvAM <- row.names(res_F0_OWvAM[res_F0_OWvAM$padj < 0.05,])

# For OA vs AM
res_F0_OAvAM <- results(dds_sub, name="treatment_OA_vs_AM", alpha=0.05)
res_F0_OAvAM <- res_F0_OAvAM[order(res_F0_OAvAM$padj),]
head(res_F0_OAvAM)

summary(res_F0_OAvAM)
res_F0_OAvAM <- res_F0_OAvAM[!is.na(res_F0_OAvAM$padj),]
degs_F0_OAvAM <- row.names(res_F0_OAvAM[res_F0_OAvAM$padj < 0.05,])

# For OWA vs AM
res_F0_OWAvAM <- results(dds_sub, name="treatment_OWA_vs_AM", alpha=0.05)
res_F0_OWAvAM <- res_F0_OWAvAM[order(res_F0_OWAvAM$padj),]
head(res_F0_OWAvAM)

summary(res_F0_OWAvAM)
res_F0_OWAvAM <- res_F0_OWAvAM[!is.na(res_F0_OWAvAM$padj),]
degs_F0_OWAvAM <- row.names(res_F0_OWAvAM[res_F0_OWAvAM$padj < 0.05,])

# Total
length(degs_F0_OAvAM)  # 602
length(degs_F0_OWvAM)  # 5517 
length(degs_F0_OWAvAM)  # 3918

# Intersections
length(intersect(degs_F0_OAvAM,degs_F0_OWvAM))  # 444
length(intersect(degs_F0_OAvAM,degs_F0_OWAvAM))  # 380
length(intersect(degs_F0_OWAvAM,degs_F0_OWvAM))  # 2743

intWA <- intersect(degs_F0_OAvAM,degs_F0_OWvAM)
length(intersect(degs_F0_OWAvAM,intWA)) # 338

# Number unique

602-444-380+338 # 116 OA
5517-444-2743+338 # 2668 OW 
3918-380-2743+338 # 1133 OWA

444-338 # 106 OA & OW
380-338 # 42 OA & OWA
2743-338 # 2405 OWA & OW


# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("OA" = 116, "OW" = 2668, "OWA" = 1133, "OA&OW" = 106, "OA&OWA" = 42, 
                "OW&OWA" = 2405, "OA&OW&OWA" = 338))


plot(fit1,  lty = 1:3, quantities = TRUE)
# lty changes the lines

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))


#cross check
1133+42+338+2405 # 3918, total OW
106+338+2405+2668 # 5517, total OWA
116+106+338+42    # 602, total OA
