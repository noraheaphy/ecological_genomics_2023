
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.
setwd("/Users/mpespeni/github/hudsonica/GO_MWU")

# Can run for each of 3 inputs, treatment groups vs AM, and can do for each GO division BP, MF, CC

# Edit these to match your data file names: 
input = read.csv("res_F0_OWvAM_LFC.csv") # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations = "trinotate_annotation_GOblastx_onlyanns_onlyGOs_justmergeGeneIDtab34.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase = "go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision = "BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.
# Can't find the super long file name that I guess I did need to transfer from FileZilla.
# Go get 'dissim0_BP_trinotate_annotation_GOblastx_onlyanns_onlyGOs_justmergeGeneIDtab34.txt'


gomwuStats(input, goDatabase, goAnnotations, goDivision,
	perlPath = "perl", # replace with full path to perl executable if it is not in your system's PATH already
	largest = 0.05,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
	smallest = 10,   # a GO category should contain at least this many genes to be considered
	clusterCutHeight = 0.25 # threshold for merging similar (gene-sharing) terms. See README for details.
#	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
#	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
#	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# OUTPUT:
# go.obo trinotate_annotation_GOblastx_onlyanns_onlyGOs_justmergeGeneIDtab34.txt res_F0_OWvAM_LFC.csv BP largest=0.1 smallest=5 cutHeight=0.25
# 
# Run parameters:
#   
#   largest GO category as fraction of all genes (largest)  : 0.1
# smallest GO category as # of genes (smallest)  : 5
# clustering threshold (clusterCutHeight) : 0.25
# 
# -----------------
#   retrieving GO hierarchy, reformatting data...
# 
# -------------
#   go_reformat:
#   Genes with GO annotations, but not listed in measure table: 72035
# 
# Terms without defined level (old ontology?..): 231
# -------------
#   -------------
#   go_nrify:
#   11061 categories, 9260 genes; size range 5-926 # max is 926 because largest=0.1 above, try making smaller
# 100 too broad
# 4832 too small
# 6129 remaining
# 
# removing redundancy:
#   
#   calculating GO term similarities based on shared genes...
# 3981 non-redundant GO categories of good size
# -------------
#   
#   Secondary clustering:
#   calculating similarities....
# Continuous measure of interest: will perform MWU test
# 853 GO terms at 10% FDR

# largest GO category as fraction of all genes (largest)  : 0.05
# smallest GO category as # of genes (smallest)  : 10
# clustering threshold (clusterCutHeight) : 0.25
# 
# -----------------
#   retrieving GO hierarchy, reformatting data...
# 
# -------------
#   go_reformat:
#   Genes with GO annotations, but not listed in measure table: 72035
# 
# Terms without defined level (old ontology?..): 231
# -------------
#   -------------
#   go_nrify:
#   11061 categories, 9260 genes; size range 10-463
# 199 too broad
# 6903 too small
# 3959 remaining

# removing redundancy:
#   
#   calculating GO term similarities based on shared genes...
# 2672 non-redundant GO categories of good size
# -------------
#   
#   Secondary clustering:
#   calculating similarities....
# Continuous measure of interest: will perform MWU test
# 700 GO terms at 10% FDR - for OW out of 2672, 26.2%

# 553 GO terms at 10% FDR - for OWA out of 2672, 20.7%

# 166 GO terms at 10% FDR - for OA out of 2672, 6.2%


# ----------- Plotting results

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
 	# absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
  absValue=1, # un-remark this if you are using log2-fold changes
 	level1=0.001, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
 	level2=0.0005, # FDR cutoff to print in regular (not italic) font.
 	level3=0.0001, # FDR cutoff to print in large bold font.
 	txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
 	treeHeight=0.5, # height of the hierarchical clustering tree
 #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
 )
 # manually rescale the plot so the tree matches the text 
 # if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  
 
png(filename = "dendroGOMWU_BP_lesssig_OAvAM_LFC.png", width = 800, height = 1400, res = 200)
# results above
dev.off()

 # text representation of results, with actual adjusted p-values
 results[[1]]


# ------- extracting representative GOs

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
plot(results[[2]],cex=0.6)
abline(h=hcut,col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
	rn=names(ct)[ct==ci]
	obs=grep("obsolete",rn)
	if(length(obs)>0) { rn=rn[-obs] }
	if (length(rn)==0) {next}
	rr=results[[1]][rn,]
	bestrr=rr[which(rr$pval==min(rr$pval)),]
	best=1
	if(nrow(bestrr)>1) {
		nns=sub(" .+","",row.names(bestrr))
		fr=c()
		for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
		best=which(fr==max(fr))
	}
	if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs

