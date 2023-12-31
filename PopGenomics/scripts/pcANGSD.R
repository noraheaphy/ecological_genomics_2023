library(ggplot2) # plotting
library(ggpubr) # plotting

setwd("~/Desktop/01_projects/PBIO_6800_ecological_genomics/ecological_genomics_2023/PopGenomics/results/") # set the path to where you saved the pcANGSD results on your laptop

## First, let's work on the genetic PCA:

COV <- as.matrix(read.table("allRS_poly.cov")) # read in the genetic covariance matrix

PCA <- eigen(COV) # extract the principal components from the COV matrix

## How much variance is explained by the first few PCs?

var <- round(PCA$values/sum(PCA$values),3)

var[1:3]

# A "screeplot" of the eigenvalues of the PCA:

barplot(var, 
        xlab = "Eigenvalues of the PCA", 
        ylab = "Proportion of variance explained")

## Bring in the bam.list file and extract the sample info:

names <- read.table("allRS_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names(pops) = c("Ind", "Pop", "Row", "Col")

## A quick and humble PCA plot:

plot(PCA$vectors[,1:2],
     col=as.factor(pops[,2]),
     xlab="PC1",ylab="PC2", 
     main="Genetic PCA")

## A more beautiful PCA plot using ggplot :)

data=as.data.frame(PCA$vectors)
data=data[,c(1:3)]
data= cbind(data, pops)

cols=c("#377eB8","#EE9B00","#0A9396","#94D2BD","#FFCB69","#005f73","#E26D5C","#AE2012", "#6d597a", "#7EA16B","#d4e09b", "gray70")

ggscatter(data, x = "V1", y = "V2",
          color = "Pop",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols), name="Source population") +
  guides(colour = guide_legend(nrow = 2))




## Next, we can look at the admixture clustering:

# import the ancestry scores (these are the .Q files)

q <- read.table("allRS_poly.admix.2.Q", sep=" ", header=F)

K=dim(q)[2] #Find the level of K modeled

## order according to population code
ord<-order(pops[,2])

# make the plot:
barplot(t(q)[,ord],
        col=cols[1:K],
        space=0,border=NA,
        xlab="Populations",ylab="Admixture proportions",
        main=paste0("Red spruce K=",K))
text(tapply(1:nrow(pops),pops[ord,2],mean),-0.05,unique(pops[ord,2]),xpd=T)
abline(v=cumsum(sapply(unique(pops[ord,2]),function(x){sum(pops[ord,2]==x)})),col=1,lwd=1.2)

####### Selection scans for red spruce ######

library(RcppCNPy) # for reading python numpy (.npy) files

setwd("~/Desktop/01_projects/PBIO_6800_ecological_genomics/ecological_genomics_2023/PopGenomics/results/")

list.files()

### read in selection statistics (these are chi^2 distributed)

s <- npyLoad("RS_mapped_poly.selection.npy")

# convert test statistic to p-value
pval <- as.data.frame(1 - pchisq(s, 1))
names(pval) = c("p_PC1", "p_PC2")

## read positions
p <- read.table("RS_mapped_poly_mafs.sites", sep = "\t", header = T, stringsAsFactors = T)
dim(p)

p_filtered = p[which(p$kept_sites==1),]
dim(p_filtered)

# How many sites got filtered out when testing for selection? Why?

## make manhattan plot
plot(-log10(pval$p_PC1),
     col = p_filtered$chromo,
     xlab = "Position",
     ylab = "-log10(p-value)",
     main = "Selection outliers: pcANGSD e = 1 (K2)")
 
# We can zoom in if there's something interesting near a position...

plot(-log10(pval$p_PC1[2e05:2.01e05]),
     col=p_filtered$chromo, 
     xlab="Position", 
     ylab="-log10(p-value)", 
     main="Selection outliers: pcANGSD e=1 (K2)")

# get the contig with the lowest p-value for selection
sel_contig <- p_filtered[which(pval==min(pval$p_PC1)), c("chromo", "position")]
sel_contig

# get all the outliers with p-values below some cutoff
cutoff = 1e-3
outliers_PC1 <- p_filtered[which(pval$p_PC1 < cutoff),c("chromo", "position")]
outliers_PC2 <- p_filtered[which(pval$p_PC2 < cutoff),c("chromo", "position")]

# how many outlier loci < the cutoff?
dim(outliers_PC1)[1]
dim(outliers_PC2)[1]

# how many unique contigs harbor outlier loci?
length(unique(outliers_PC1$chromo))

# export outlier contigs
write.table(outliers_PC1,
            "allRS_poly_outliers_PC1.txt", 
            sep = ":",
            quote = F,
            row.names = F,
            col.names = F)

COV <- as.matrix(read.table("allRS_poly.cov"))

PCA <- eigen(COV)

data = as.data.frame(PCA$vectors)
data = data[,c(1:2)] # the second number here is the number of PC axes you want to keep

write.table(data,
            "allRS_poly_genPC1_2.txt",
            sep="\t",
            quote = F,
            row.names = F,
            col.names = F)
