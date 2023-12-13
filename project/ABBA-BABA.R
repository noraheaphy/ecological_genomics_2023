library(tidyverse)

# read in ABBA-BABA results
genome_wide <-  read.csv("result.Observed.txt", header = T, sep = "\t", na.strings = "")
combined <- read.csv("combined_file_v2.txt", header = F, sep = "\t", na.strings = "")

# add in column names
column_names <- c("contig", "D", "JK_D", "V_JK_D", "Z", "pvalue", "nABBA", "nBABA", "nBlocks", 
                  "H1", "H2", "H3", "H4")
colnames(combined) <- column_names

# fix contig names in first column
combined$contig <- str_replace(combined$contig, "_row5.txt", "")

# correct for multiple testing
combined$p_adjust <- p.adjust(combined$pvalue, method = "BH")


# pull out significant ones
signif <- combined %>% filter(p_adjust <= 0.0005)

# save to file
write.csv(signif, "signif_contigs_abba_v2.csv", row.names = F)
write.csv(signif$contig, "signif_contigs_name_only.txt", row.names = F)

## pull in PC loadings significant genes (from original mapping)
pc_signif <- read.csv("RSBS_PC1_outlier_genes.txt", header = F)
pc_signif <- pc_signif$V1

# pull in genes significant from local PCA
local_signif <- read.csv("local_PCA_geneIDs.txt", header = F)
local_signif <- local_signif$V1

# pull in gene_IDs from ABBA-BABA signif
abba_genes <- read.csv("gene_IDs.txt", header = F)

overlap_pc_abba <- abba_genes %>% filter(V1 %in% pc_signif)
overlap_abba_local <- abba_genes %>% filter(V1 %in% local_signif)

overlap_pc_local <- pc_signif %>% filter(V1 %in% local_signif)
