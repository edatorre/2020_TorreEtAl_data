##The goal of this script is to obtain the differentially expressed genes from RNA seq data from PreResistant cells


#Load datasets
#This is the raw counts file
t01_allData = read.table("01_rawData/190408_RNAseq_otherSamples/NGFRandEGFRhi/rawCounts/meltedData.tsv", header = T)

#Convert tall table into matrix
allCounts_rawCountsMatrix = t01_allData %>%
  spread(sampleID, counts)

#Clean up table
allCounts_rawCountsMatrix = allCounts_rawCountsMatrix %>%
  select(-experiment)

allCounts_rawCountsMatrix = allCounts_rawCountsMatrix[!grepl("_", allCounts_rawCountsMatrix$gene_id),]

#Add row names
rownames(allCounts_rawCountsMatrix) = allCounts_rawCountsMatrix$gene_id
allCounts_rawCountsMatrix$gene_id = NULL

#Rearrange columns
allCounts_rawCountsMatrix = allCounts_rawCountsMatrix %>%
  select(`PreR-001-BE-neg`, `PreR-003-ET-neg`, `PreR-005-YG-neg`, `PreR-002-BE-DP`, `PreR-004-ET-DP`, `PreR-006-YG-DP`)

#Set conditions and format them 
condition = factor(c("wt", "wt", "wt", "dp", "dp", "dp"), levels=c("wt", "dp")) ###This is what determines the order of the comparison.
coldata = data.frame(condition)
coldata$type = factor("paired-end")
rownames(coldata) = colnames(allCounts_rawCountsMatrix)

#Generate DEseq data set.
dds <- DESeqDataSetFromMatrix(countData = allCounts_rawCountsMatrix,
                              colData = coldata,
                              design = ~ condition)
#Run DEseq
dds <- DESeq(dds)

#Extract Results
res <- results(dds)
res = data.frame(res)
res$id = rownames(res)

#filter to keep differentially expressed genes
res_filtered = res %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>%
  select(id, log2FoldChange)


#Save the tables
write.table(res_filtered, file = "03_extractedData/190408_RNAseq_otherSamples/NGFRandEGFRhi/DEseq/NGFRhiEGFRhi_differentiallyExpressedGenes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(res, file = "03_extractedData/190408_RNAseq_otherSamples/NGFRandEGFRhi/DEseq/NGFRhiEGFRhi_allDifferentialExpression.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



