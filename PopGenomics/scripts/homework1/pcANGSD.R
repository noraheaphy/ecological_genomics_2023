library(ggplot2)
library(ggpubr) 

setwd("~/Desktop/01_projects/PBIO_6800_ecological_genomics/ecological_genomics_2023/PopGenomics/results/homework1") # set the path to where you saved the pcANGSD results on your laptop

## Bring in the bam.list file and extract the sample info

names <- read.table("allRS_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names(pops) = c("Ind", "Pop", "Row", "Col")

# set up colors ahead of time
cols = c("#377eB8","#EE9B00","#0A9396","#94D2BD","#FFCB69","#005f73","#E26D5C","#AE2012", "#6d597a", "#7EA16B","#d4e09b", "gray70")

#### K = 2 -----------------------------------------

## First, let's work on the genetic PCA:

COV2 <- as.matrix(read.table("allRS_poly2.cov")) # read in the genetic covariance matrix

PCA2 <- eigen(COV2) # extract the principal components from the COV matrix

## How much variance is explained by the first few PCs?

var2 <- round(PCA2$values/sum(PCA2$values), 3)

var2[1:3]

# A "screeplot" of the eigenvalues of the PCA:

barplot(var2, 
        main = "K=2",
        xlab = "Eigenvalues of the PCA", 
        ylab = "Proportion of variance explained",
        cex.lab = 1.5, cex.axis = 1.2, cex.main = 2)

## A quick and humble PCA plot:

plot(PCA2$vectors[,1:2],
     col = as.factor(pops[,2]),
     xlab = "PC1", ylab = "PC2", 
     main = "Genetic PCA")

## A more beautiful PCA plot using ggplot :)

data2 = as.data.frame(PCA2$vectors)
data2 = data2[,c(1:3)]
data2 = cbind(data2, pops)

ggscatter(data2, x = "V1", y = "V2",
          color = "Pop",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13) +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(.7)), 
        axis.text = element_text(size = 13), 
        legend.position = "left") +
  labs(x = paste0("PC1: (",var2[1]*100,"%)"), y = paste0("PC2: (",var2[2]*100,"%)")) +
  scale_color_manual(values = c(cols), name = "Source population") +
  guides(colour = guide_legend(nrow = 6))

## Next, we can look at the admixture clustering:

# import the ancestry scores (these are the .Q files)

q2 <- read.table("allRS_poly2.admix.2.Q", sep = " ", header = F)

K2 = dim(q2)[2] #Find the level of K modeled

## order according to population code
ord2 <- order(pops[,2])

# make the plot:
barplot(t(q2)[,ord2],
        col = cols[1:K2],
        space = 0, border = NA,
        xlab = "Populations", ylab = "Admixture proportions",
        main = paste0("Red spruce K=", K2))
text(tapply(1:nrow(pops), pops[ord2,2], mean), -0.05, unique(pops[ord2,2]), xpd = T)
abline(v = cumsum(sapply(unique(pops[ord2,2]), function(x){sum(pops[ord2,2]==x)})), col = 1, lwd = 1.2)

#### K = 3 -----------------------------------------

## First, let's work on the genetic PCA:

COV3 <- as.matrix(read.table("allRS_poly3.cov")) # read in the genetic covariance matrix

PCA3 <- eigen(COV3) # extract the principal components from the COV matrix

## How much variance is explained by the first few PCs?

var3 <- round(PCA3$values/sum(PCA3$values), 3)

var3[1:3]

# A "screeplot" of the eigenvalues of the PCA:

barplot(var3, 
        main = "K=3",
        xlab = "Eigenvalues of the PCA", 
        ylab = "Proportion of variance explained",
        cex.lab = 1.5, cex.axis = 1.2, cex.main = 2)

## A quick and humble PCA plot:

plot(PCA3$vectors[,1:2],
     col = as.factor(pops[,2]),
     xlab = "PC1", ylab = "PC2", 
     main = "Genetic PCA")

## A more beautiful PCA plot using ggplot :)

data3 = as.data.frame(PCA3$vectors)
data3 = data3[,c(1:3)]
data3 = cbind(data3, pops)

ggscatter(data3, x = "V1", y = "V2",
          color = "Pop",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13) +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(.7)), 
        axis.text = element_text(size = 13), 
        legend.position = "left") +
  labs(x = paste0("PC1: (",var3[1]*100,"%)"), y = paste0("PC2: (",var3[2]*100,"%)")) +
  scale_color_manual(values = c(cols), name = "Source population") +
  guides(colour = guide_legend(nrow = 6))

## Next, we can look at the admixture clustering:

# import the ancestry scores (these are the .Q files)

q3 <- read.table("allRS_poly3.admix.3.Q", sep = " ", header = F)

K3 = dim(q3)[2] #Find the level of K modeled

## order according to population code
ord3 <- order(pops[,2])

# make the plot:
barplot(t(q3)[,ord3],
        col = cols[1:K3],
        space = 0, border = NA,
        xlab = "Populations", ylab = "Admixture proportions",
        main = paste0("Red spruce K=", K3))
text(tapply(1:nrow(pops), pops[ord3,2], mean), -0.05, unique(pops[ord3,2]), xpd = T)
abline(v = cumsum(sapply(unique(pops[ord3,2]), function(x){sum(pops[ord3,2]==x)})), col = 1, lwd = 1.2)

#### K = 4 -----------------------------------------

## First, let's work on the genetic PCA:

COV4 <- as.matrix(read.table("allRS_poly4.cov")) # read in the genetic covariance matrix

PCA4 <- eigen(COV4) # extract the principal components from the COV matrix

## How much variance is explained by the first few PCs?

var4 <- round(PCA4$values/sum(PCA4$values), 3)

var4[1:3]

# A "screeplot" of the eigenvalues of the PCA:

barplot(var4, 
        main = "K=4",
        xlab = "Eigenvalues of the PCA", 
        ylab = "Proportion of variance explained",
        cex.lab = 1.5, cex.axis = 1.2, cex.main = 2)

## A quick and humble PCA plot:

plot(PCA4$vectors[,1:2],
     col = as.factor(pops[,2]),
     xlab = "PC1", ylab = "PC2", 
     main = "Genetic PCA")

## A more beautiful PCA plot using ggplot :)

data4 = as.data.frame(PCA4$vectors)
data4 = data4[,c(1:3)]
data4 = cbind(data4, pops)

ggscatter(data4, x = "V1", y = "V2",
          color = "Pop",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13) +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(.7)), 
        axis.text = element_text(size = 13), 
        legend.position = "left") +
  labs(x = paste0("PC1: (",var4[1]*100,"%)"), y = paste0("PC2: (",var4[2]*100,"%)")) +
  scale_color_manual(values = c(cols), name = "Source population") +
  guides(colour = guide_legend(nrow = 6))

## Next, we can look at the admixture clustering:

# import the ancestry scores (these are the .Q files)

q4 <- read.table("allRS_poly4.admix.4.Q", sep = " ", header = F)

K4 = dim(q4)[2] #Find the level of K modeled

## order according to population code
ord4 <- order(pops[,2])

# make the plot:
barplot(t(q4)[,ord4],
        col = cols[1:K4],
        space = 0, border = NA,
        xlab = "Populations", ylab = "Admixture proportions",
        main = paste0("Red spruce K=", K4))
text(tapply(1:nrow(pops), pops[ord4,2], mean), -0.05, unique(pops[ord4,2]), xpd = T)
abline(v = cumsum(sapply(unique(pops[ord4,2]), function(x){sum(pops[ord4,2]==x)})), col = 1, lwd = 1.2)


#### selection scans K = 3 --------------------

library(RcppCNPy) # for reading python numpy (.npy) files

setwd("~/Desktop/01_projects/PBIO_6800_ecological_genomics/ecological_genomics_2023/PopGenomics/results/homework1")

list.files()

### read in selection statistics (these are chi^2 distributed)

s <- npyLoad("allRS_poly.selection.npy")

# convert test statistic to p-value
pval <- as.data.frame(1 - pchisq(s, 1))
names(pval) = c("p_PC1", "p_PC2")

## read positions
p <- read.table("allRS_poly_mafs.sites", sep = "\t", header = T, stringsAsFactors = T)
dim(p)

p_filtered = p[which(p$X1==1),]
dim(p_filtered)

# How many sites got filtered out when testing for selection? Why?

## make manhattan plot
plot(-log10(pval$p_PC1),
     col = p_filtered$chromo,
     xlab = "Position",
     ylab = "-log10(p-value)",
     main = "Selection outliers: pcANGSD e = 2 (K3)")

# We can zoom in if there's something interesting near a position...

plot(-log10(pval$p_PC1[2e05:2.01e05]),
     col = p_filtered$chromo, 
     xlab = "Position", 
     ylab = "-log10(p-value)", 
     main = "Selection outliers: pcANGSD e = 2 (K3)")

# get the contig with the lowest p-value for selection
sel_contig <- p_filtered[which(pval==min(pval$p_PC1)), c("chromo", "position")]
sel_contig

# get all the outliers with p-values below some cutoff
cutoff3 = 1e-3
outliers_PC1_cutoff3 <- p_filtered[which(pval$p_PC1 < cutoff3), c("chromo", "position")]
outliers_PC2_cutoff3 <- p_filtered[which(pval$p_PC2 < cutoff3), c("chromo", "position")]

cutoff4 = 1e-4
outliers_PC1_cutoff4 <- p_filtered[which(pval$p_PC1 < cutoff4), c("chromo", "position")]
outliers_PC2_cutoff4 <- p_filtered[which(pval$p_PC2 < cutoff4), c("chromo", "position")]

# how many outlier loci < the cutoff?
dim(outliers_PC1)[1]
dim(outliers_PC2)[1]

# how many unique contigs harbor outlier loci?
length(unique(outliers_PC1$chromo))

# export outlier contigs
write.table(outliers_PC1_cutoff3,
            "allRS_poly_outliers_PC1_cutoff3.txt", 
            sep = ":",
            quote = F,
            row.names = F,
            col.names = F)

write.table(outliers_PC2_cutoff3,
            "allRS_poly_outliers_PC2_cutoff3.txt", 
            sep = ":",
            quote = F,
            row.names = F,
            col.names = F)

write.table(outliers_PC1_cutoff4,
            "allRS_poly_outliers_PC1_cutoff4.txt", 
            sep = ":",
            quote = F,
            row.names = F,
            col.names = F)

write.table(outliers_PC2_cutoff4,
            "allRS_poly_outliers_PC2_cutoff4.txt", 
            sep = ":",
            quote = F,
            row.names = F,
            col.names = F)

COV <- as.matrix(read.table("allRS_poly3.cov"))

PCA <- eigen(COV)

data = as.data.frame(PCA$vectors)
data = data[,c(1:2)] # the second number here is the number of PC axes you want to keep

write.table(data,
            "allRS_poly_genPC1_2.txt",
            sep="\t",
            quote = F,
            row.names = F,
            col.names = F)

