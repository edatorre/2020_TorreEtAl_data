# #The goal of this script is to run enrichment analysi on the leading edges of different clusters

listOfFiles = list.files("03_extractedData/190321_KO-RNAseq/GSEA/leadingEdges/")

#Load transcriptome and keep only genes expressed in at least one sample.
expressedTranscriptome = read.table("01_rawData/190321_KO-RNAseq/normalizedCounts/rpmCounts_allRuns_matrix.tsv", header = T)
expressedTranscriptome[is.na(expressedTranscriptome)] = 0
expressedTranscriptome = expressedTranscriptome[rowSums(expressedTranscriptome[3:ncol(expressedTranscriptome)]) > 0, ]

#generate background and convert names
background = as.character(expressedTranscriptome$GeneSymbol)
background = bitr(background, fromType="SYMBOL", toType=("ENTREZID"), OrgDb="org.Hs.eg.db")

for(currentFile in listOfFiles) {
  
  # #For testing
  # currentFile = listOfFiles[2]
  
  #Remove ".txt" to use base of name for file names
  currentFile_baseName = strsplit(currentFile, ".txt")
  
  #Load leading edge genes
  leadingEdge = read.table(paste("03_extractedData/190321_KO-RNAseq/GSEA/leadingEdges/", currentFile, sep = ""), sep = "\t")
  
  #convert names of the leading edge
  leadingEdge = bitr(leadingEdge$V1, fromType="SYMBOL", toType=("ENTREZID"), OrgDb="org.Hs.eg.db")
  
  #compute enrichment analysis - upregulated genes
  ego1 <- enrichGO(gene          = leadingEdge$ENTREZID,
                   universe      = background$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 0.25,
                   readable      = TRUE)
  ego1.df = data.frame(ego1)
  rownames(ego1.df) = NULL
  
  newFileName_GO = paste("03_extractedData/190321_KO-RNAseq/GSEA/enrichmentAnalysisFromClusters/", currentFile_baseName[[1]][1], "_GOenrichment.txt", sep = "")
  write.table(ego1.df, file = newFileName_GO, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  kk <- enrichKEGG(gene         = leadingEdge$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.25)
  kk.df = data.frame(kk)
  
  newFileName_KEGG = paste("03_extractedData/190321_KO-RNAseq/GSEA/enrichmentAnalysisFromClusters/", currentFile_baseName[[1]][1], "_KEGGenrichment.txt", sep = "")
  write.table(kk.df, file = newFileName_KEGG, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  
}
