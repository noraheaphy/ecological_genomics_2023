## wd is "/Users/nheaphy/Desktop/01_projects/PBIO_6800_ecological_genomics/ecological_genomics_2023/Transcriptomics/GOMWU"

## Import the libraries that we're likely to need in this session
library(DESeq2)

# Try with new counts table from filtered transcriptome assembly
countsTable <- read.table("../myresults/salmon.isoform.counts.matrix.filteredAssembly", header = TRUE, row.names = 1)

head(countsTable)
dim(countsTable)
#[1] 130580     38 - genes
# [1] 349516     38 - isoforms
# [1] 67916    38 - filtered assembly

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample discription table
conds <- read.delim("../myresults/ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

#################### MODEL NUMBER 2 - subset to focus on effect of treatment for each generation

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ treatment)

dim(dds)
# [1] 130580     38

# Filter 
dds <- dds[rowSums(counts(dds) >= 30) >= 28,]
nrow(dds) 

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)

# Perform DESeq2 analysis on the subset
dds_F0 <- DESeq(dds_F0)

resultsNames(dds_F0)
# [1] "Intercept"           "treatment_OA_vs_AM"  "treatment_OW_vs_AM"  "treatment_OWA_vs_AM"

res_F0_OWvAM <- results(dds_F0, name = "treatment_OW_vs_AM", alpha = 0.05)

res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM) 

summary(res_F0_OWvAM)

res_F0_OWAvAM <- results(dds_F0, name = "treatment_OWA_vs_AM", alpha = 0.05)

res_F0_OWAvAM <- res_F0_OWAvAM[order(res_F0_OWAvAM$padj),]
head(res_F0_OWAvAM) 

summary(res_F0_OWAvAM)


res_F0_OAvAM <- results(dds_F0, name = "treatment_OA_vs_AM", alpha = 0.05)

res_F0_OAvAM <- res_F0_OAvAM[order(res_F0_OAvAM$padj),]
head(res_F0_OAvAM)

summary(res_F0_OAvAM)


################## Save all the results as csv to go into GOMWU

library(tidyr)

# Make the rownames a separate column called transcriptID and make it all a dataframe
res_F0_OWvAM_df <- data.frame(transcriptID = rownames(res_F0_OWvAM), res_F0_OWvAM)

# Split the "transcriptID" column by double colons and create new columns of the parts
res_F0_OWvAM_df <- separate(res_F0_OWvAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 

# Create a new column by concatenating "part1" and "part2" with double colons in between
res_F0_OWvAM_df$transcriptID_trim <- paste(res_F0_OWvAM_df$part1, res_F0_OWvAM_df$part2, sep = "::")

# Optional: Remove the "part1" and "part2" columns from the dataframe
res_F0_OWvAM_df <- res_F0_OWvAM_df[, !(names(res_F0_OWvAM_df) %in% c("part1", "part2", "part3", "rest"))]

write.table(res_F0_OWvAM_df, file = "res_F0_OWvAM.txt", sep = "\t", row.names = F)   # saves the full original for the records

# Select the two columns we want to save for the GOMWU analysis
selected_columns_OW <- res_F0_OWvAM_df[c("transcriptID_trim", "log2FoldChange")]

# Save the selected columns as a CSV file
write.csv(selected_columns_OW, file = "res_F0_OWvAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU


############ Now for OWA

# Make the rownames a separate column called transcriptID and make it all a dataframe
res_F0_OWAvAM_df <- data.frame(transcriptID = rownames(res_F0_OWAvAM), res_F0_OWAvAM)

# Split the "transcriptID" column by double colons and create new columns of the parts
res_F0_OWAvAM_df <- separate(res_F0_OWAvAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 

# Create a new column by concatenating "part1" and "part2" with double colons in between
res_F0_OWAvAM_df$transcriptID_trim <- paste(res_F0_OWAvAM_df$part1, res_F0_OWAvAM_df$part2, sep = "::")

# Optional: Remove the "part1" and "part2" columns from the dataframe
res_F0_OWAvAM_df <- res_F0_OWAvAM_df[, !(names(res_F0_OWAvAM_df) %in% c("part1", "part2", "part3", "rest"))]
write.table(res_F0_OWAvAM_df, file = "res_F0_OWAvAM.txt", sep = "\t", row.names = F)   # saves the full original for the records

# Select the two columns we want to save for the GOMWU analysis
selected_columns_OWA <- res_F0_OWAvAM_df[c("transcriptID_trim", "log2FoldChange")]

# Save the selected columns as a CSV file
write.csv(selected_columns_OWA, file = "res_F0_OWAvAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU



############ Now for OA

# Make the rownames a separate column called transcriptID and make it all a dataframe
res_F0_OAvAM_df <- data.frame(transcriptID = rownames(res_F0_OAvAM), res_F0_OAvAM)

# Split the "transcriptID" column by double colons and create new columns of the parts
res_F0_OAvAM_df <- separate(res_F0_OAvAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 

# Create a new column by concatenating "part1" and "part2" with double colons in between
res_F0_OAvAM_df$transcriptID_trim <- paste(res_F0_OAvAM_df$part1, res_F0_OAvAM_df$part2, sep = "::")

# Optional: Remove the "part1" and "part2" columns from the dataframe
res_F0_OAvAM_df <- res_F0_OAvAM_df[, !(names(res_F0_OAvAM_df) %in% c("part1", "part2", "part3", "rest"))]
write.table(res_F0_OAvAM_df, file = "res_F0_OAvAM.txt", sep = "\t", row.names = F)   # saves the full original for the records

# Select the two columns we want to save for the GOMWU analysis
selected_columns_OA <- res_F0_OAvAM_df[c("transcriptID_trim", "log2FoldChange")]

# Save the selected columns as a CSV file
write.csv(selected_columns_OA, file = "res_F0_OAvAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU
