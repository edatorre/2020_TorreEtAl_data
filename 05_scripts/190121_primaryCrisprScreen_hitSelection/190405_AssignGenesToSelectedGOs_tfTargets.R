# #Load libraries that will be used in the scripts.
# library(tidyverse)
# 
# #clear workspace
# rm(list=ls())

#Go terms
GOs = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/GOenrichment_TFLib_manuallyCurated.txt", header = T, sep = "\t")

GOs = GOs %>%
  filter(Keep == "yes")

#Screen tagrets
screenTargets = read.table("02_metadata/ScreenTargets_metadata.txt", header = T)

tfTargets = screenTargets %>%
  filter(library == "TF") %>%
  dplyr::select(gene) %>%
  distinct()

nGO = nrow(GOs)
outputFile = data.frame(gene = as.character(""), ID = as.character(""), Description = as.character(""))
outputFile = outputFile[-1,]

#Now loop through the GO terms and fill in the summary matrix
for(i in 1:nGO) {
  
  #for testing
  #i = 8
  
  mainGO = GOs$ID[i]
  mainGO_description = GOs$Description[i]
  message = paste(i, "out of", nGO, "-", mainGO, sep = " ")
  print(message)
  flush.console()
  
  #rownames(summaryMatrix)[i] = as.character(mainGO_description)
  
  #Make table with genes in GO
  genesInGO_string = as.character(GOs$geneID[i])
  genesInGO <- data.frame(strsplit(genesInGO_string, "/"))
  colnames(genesInGO) = "GeneSymbol"
  
  #Genes to add to output file
  genesForOutput = tfTargets %>%
    filter(gene %in% genesInGO$GeneSymbol)
  
  if(nrow(genesForOutput) > 0) {
    
  genesForOutput$ID = mainGO
  genesForOutput$Description = mainGO_description
  
  outputFile = bind_rows(outputFile, genesForOutput)
  
  #Eliminate those genes from the tfTargets file
  tmp = tfTargets %>%
    filter(!gene %in% genesForOutput$gene)
  
  #replace tfTargets with tmp
  tfTargets = tmp
   
  }
  
}
  
tfTargets$Description = "other"

outputFile = bind_rows(outputFile, tfTargets)

setwd("03_extractedData/190121_primaryCrisprScreen_hitSelection/")
write.table(outputFile, "TFTargets_assignedToGOs.txt", sep = "\t", col.names = T, row.names = FALSE, quote = FALSE)
