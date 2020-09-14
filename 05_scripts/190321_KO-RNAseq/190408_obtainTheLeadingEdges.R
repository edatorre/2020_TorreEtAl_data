# #The goal of this script is to extract the leading edges of the gene sets in a given cluster

#Load cluster contents
clusterContents_geneSets = read.table("03_extractedData/190321_KO-RNAseq/GSEA/clusterContents/c5.bp.v6.2_clusterContents_geneSets.txt", header = T)
clusterContents_KOs = read.table("03_extractedData/190321_KO-RNAseq/GSEA/clusterContents/c5.bp.v6.2_clusterContents_KOs.txt", header = T)

#Add LtoR order of gene set clusters
tmp = clusterContents_geneSets %>%
  dplyr::select(col_dividedTree) %>%
  distinct()
tmp = tmp %>%
  mutate(LtoR = 1:nrow(tmp))
clusterContents_geneSets = left_join(clusterContents_geneSets, tmp, by = "col_dividedTree")

#Add TtoB order of KO clusters
tmp = clusterContents_KOs %>%
  dplyr::select(row_dividedTree) %>%
  distinct()
tmp = tmp %>%
  mutate(TtoB = 1:nrow(tmp))
clusterContents_KOs = left_join(clusterContents_KOs, tmp, by = "row_dividedTree")

#Obtain a list of gene sets
listOfGeneSets = clusterContents_geneSets %>%
  filter(LtoR == selectedGeneSetCluster)
listOfGeneSets = as.character(listOfGeneSets$GS)

#Obtain a list of KOs in a specific cluster
listOfKO = clusterContents_KOs %>%
  filter(TtoB == selectedKOcluster) 
listOfKO = as.character(listOfKO$KO)

#Make empty data.frame to store genes
leadingEdges = data.frame(gene = as.character(), KO = as.character())

#Loop through each of the KOs
for(ko in listOfKO) {
  
  # #For testing
  # ko = listOfKO[12]
  
  #Navigate to the directory with the GSEA output for the KO
  pathToGSEAfiles = paste("03_extractedData/190321_KO-RNAseq/GSEA/GSEAoutput/c5.bp.v6.2/", ko, "/", sep = "")
  
  #obtain list of files
  listOfFiles = list.files(pathToGSEAfiles)
  listOfFiles = as.character(listOfFiles)
  
  #Remove all PDFs
  listOfFiles = listOfFiles[grepl("*report*", listOfFiles)]
  listOfFiles = listOfFiles[!grepl("*SUMMARY.RESULTS.REPORT*", listOfFiles)] #maybe not necessary
  
  #Loop thorugh the gene sets in the cluster
  for (geneSet in listOfGeneSets) {
     
    # #for testing
    # geneSet = listOfGeneSets[6]
    
    #select file form list
    currentGeneSetFiles_Name = listOfFiles[grepl(geneSet, listOfFiles)]
    
    #Check if the file exists. 
    if(!is_empty(currentGeneSetFiles_Name)) {
      
      for (i in currentGeneSetFiles_Name) {
        
        #upload file 
        currentGeneSetFile = read.table(paste(pathToGSEAfiles, i, sep = "")) #######im here#####
        
        #Select leading edge genes
        leadingEdgeGenes = currentGeneSetFile %>%
          filter(V8 == "YES") %>%
          dplyr::select(V2) %>%
          mutate(KO = ko)
        
        colnames(leadingEdgeGenes) = c("gene", "KO")
        
        leadingEdges = bind_rows(leadingEdges, leadingEdgeGenes)
        
        
      }
      
      
    }
    
    
  }
  

}

#Clean up leading edge file
leadingEdges = leadingEdges %>%
  distinct()

leadingEdges_genesOnly = leadingEdges %>%
  dplyr::select(gene) %>%
  distinct()

#Save leading edge File
write.table(leadingEdges_genesOnly, file = paste("03_extractedData/190321_KO-RNAseq/GSEA/leadingEdges/leadingEdges_GeneSetCluster", selectedGeneSetCluster, "_KOcluster", selectedKOcluster, ".txt", sep = ""),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
