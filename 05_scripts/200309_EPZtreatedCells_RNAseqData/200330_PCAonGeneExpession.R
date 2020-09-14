# #The goal of this script is to do a PCA anlysis on a select number of genes
# 


#Load name conversion key
nameConversionKey = read.table("02_metadata/190321_KO-RNAseq/hg19_nameConvertionFile.txt", header = T)


#Load metadata
metadata = read.table("02_metadata/200309_EPZtreatedCells_RNAseqData/200312_PCAmetadata.txt", header = T)
metadata$sample = NULL
colnames(metadata)[2] = "sampleName"

#Load enrichment score matrix
ESmatrix = read.table("01_rawData/200309_EPZtreatedCells_RNAseqData/normalizedConts/rpmCounts_matrix.tsv", header = TRUE)

#Remove any rows that do not correcpong to genes
ESmatrix = ESmatrix[!grepl("__", ESmatrix$gene_id),]

#Eliminate gene symbol column and add rownames
rownames(ESmatrix) = ESmatrix$gene_id
ESmatrix$gene_id = NULL

ESmatrix = ESmatrix[rowSums(abs(ESmatrix[1:ncol(ESmatrix)])) > 18 ,]

#Make copy of ESmatrix to read later
ESmatrix_readable = ESmatrix
ESmatrix_readable$gene_id = rownames(ESmatrix_readable)
ESmatrix_readable = left_join(ESmatrix_readable, nameConversionKey, by = "gene_id")
rownames(ESmatrix_readable) = ESmatrix_readable$GeneSymbol_gene_id

#Transpose 
ESmatrix = t(ESmatrix)
# enrichmentScoreMatrix_newSamples = t(enrichmentScoreMatrix_newSamples)

#Do PCA
PCAresults = prcomp(ESmatrix, scale. = TRUE)

#Extract x
t03_XfromPCA = as.data.frame(PCAresults$x)

#remove row names and add field for the names
t03_XfromPCA$sampleName = rownames(t03_XfromPCA)
rownames(t03_XfromPCA) = NULL

#Combine metaData with PCA results
t03_XfromPCA = left_join(t03_XfromPCA, metadata, by = c("sampleName"))

#Summary on PCA results
t04 = summary(PCAresults)

#save info regarding proportion of variance accounted for by each PC
t04_importance = data.frame(t04$importance)

#Obtain loadings
loadings = data.frame(PCAresults$rotation)
loadings$gene_id = rownames(loadings)
loadings = left_join(loadings, nameConversionKey, by = "gene_id")
rownames(loadings) = loadings$GeneSymbol_gene_id

#save loadings
write.table(loadings, file = "03_extractedData/200309_EPZtreatedCells_RNAseqData/PCA/200319_loadings.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# #Plot PCA results: PC1 vs PC3 by sample type
# tpmPCAplot_test <- ggplot(t03_XfromPCA, aes(PC1, PC3)) +
#   geom_point(aes(color = sampleType)) +
#   geom_text_repel(data = t03_XfromPCA, aes(PC1, PC3, label = sampleName, color = sampleType)) +
#   theme_classic()
# plot(tpmPCAplot_test)
# ggsave("04_plots/200309_EPZtreatedCells_RNAseqData/PCA_GeneExpMatrix_PC1andPC3_byType.PDF", width = 10, height = 7)

#Plot PCA results: PC1 and PC3 by phenotype
tpmPCAplot_test <- ggplot(t03_XfromPCA, aes(PC1, PC3)) +
  geom_point(aes(color = phenotype)) +
  geom_text_repel(data = t03_XfromPCA, aes(PC1, PC3, label = sampleName, color = phenotype)) +
  theme_classic()
plot(tpmPCAplot_test)
ggsave("04_plots/200309_EPZtreatedCells_RNAseqData/PCA_GeneExpMatrix_PC1andPC3_byPhenotype.PDF", width = 10, height = 7)

#Plot PCA results: PC1 vs PC3 by pre-treatment regimen
tpmPCAplot_test <- ggplot(t03_XfromPCA, aes(PC1, PC3)) +
  geom_point(aes(color = preTregimen)) +
  geom_text_repel(data = t03_XfromPCA, aes(PC1, PC3, label = sampleName, color = preTregimen)) +
  theme_classic()
plot(tpmPCAplot_test)
ggsave("04_plots/200309_EPZtreatedCells_RNAseqData/PCA_GeneExpMatrix_PC1andPC3_byPreTregimen.PDF", width = 10, height = 7)

#Plot PCA results: PC1 vs PC2 by phenotype
tpmPCAplot_test <- ggplot(t03_XfromPCA, aes(PC1, PC2)) +
  geom_point(aes(color = phenotype)) +
  geom_text_repel(data = t03_XfromPCA, aes(PC1, PC2, label = sampleName, color = phenotype)) +
  theme_classic()
plot(tpmPCAplot_test)
ggsave("04_plots/200309_EPZtreatedCells_RNAseqData/PCA_GeneExpMatrix_PC1andPC2_byPhenotype.PDF", width = 10, height = 7)

#Plot PCA results: PC1 vs PC2 by pre-treatment regimen
tpmPCAplot_test <- ggplot(t03_XfromPCA, aes(PC1, PC2)) +
  geom_point(aes(color = preTregimen)) +
  geom_text_repel(data = t03_XfromPCA, aes(PC1, PC2, label = sampleName, color = preTregimen)) +
  theme_classic()
plot(tpmPCAplot_test)
ggsave("04_plots/200309_EPZtreatedCells_RNAseqData/PCA_GeneExpMatrix_PC1andPC2_byPreTregimen.PDF", width = 10, height = 7)




















