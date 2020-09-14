#The goal of this script is to perform GSEA in each of the KOs using the NGFR clusters as the gene sets

#Select controls to eliminate from analysis
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
negativeControls = c("neg", "WT")

#Load gene expresion data
geneExpData = read.table("01_rawData/190321_KO-RNAseq/normalizedCounts/rpmCounts_allRuns_matrix.tsv", header = T)

#Clean up the table
geneExpData[is.na(geneExpData)] = 0 #convert NAs to 0
geneExpData = geneExpData[rowSums(geneExpData[3:ncol(geneExpData)]) > 0, ] #eliminate genes not expressed at all.

# This tells me which genes remain duplicated, which I remove.
tmp = geneExpData %>%
  select(GeneSymbol) %>%
  mutate(dup = duplicated(geneExpData$GeneSymbol)) %>%
  filter(dup == T)
tmp = tmp %>%
  distinct()

removeDuplicatedGene = tmp$GeneSymbol

#Remove the duplicated gene
geneExpData = geneExpData %>%
  filter(!GeneSymbol %in% removeDuplicatedGene) %>%
  distinct()

#Now convert the gene Symbol into the row names and eliminate unnecessary fields
row.names(geneExpData) = geneExpData$GeneSymbol #add rownames to prevent losing them
geneExpData$gene_id = NULL #eliminate column
geneExpData$GeneSymbol = NULL #eliminate column

# #Load gene name conversion files
# nameConversionFile = read.table("../../metadata/hg19_nameConvertionFile.txt", header = T)

#create the first two columns of the output file
geneNmes.df = data.frame(rownames(geneExpData))
colnames(geneNmes.df) = "gene_id"

geneNmes.df = geneNmes.df %>%
  mutate(DESCRIPTION = "na")

#create the table of controls
control.geneExp = geneExpData %>%
  select(contains("neg"))

control.geneExp = control.geneExp[, colSums(control.geneExp) > 0]
control.geneExp$gene_id = rownames(control.geneExp)

#Load list of KOs
listOfKO = read.table("02_metadata/TargetedScreen_targets.txt", header = T)
listOfKO = listOfKO %>%
  dplyr::select(target) %>%
  distinct()
listOfKO = listOfKO %>%
  filter(!target %in% positiveControls) %>%
  filter(!target %in% negativeControls)

#fix nomenclature
listOfKO = gsub("-", ".", listOfKO$target) #substitute _ for . (that is the separator in the gene expression dataset)

# remove = c("POLR2A", "TP73", "WT") #these two KOs are making the whole thing crash. Im removing them for now. 
# listOfKO = listOfKO[!listOfKO %in% remove]

for (i in listOfKO) {
  
  print(i)
  flush.console()
  
  #Select KO of interest
  KO = i
  #KO = "KDM1A"
  
  #Select columns of interest
  ko.geneExp = geneExpData %>%
    select(contains(KO))
  
  #Keep only columns that have content
  ko.geneExp = ko.geneExp[, colSums(ko.geneExp) > 0]
  ko.geneExp = data.frame(ko.geneExp)
  
  #Determines the number of biological replicates in a given KO
  replicatesPerKO = ncol(ko.geneExp)  #Im eliminating samples where I dont hve a replicate. 
  
  #Make GCT file and save it. BUT, keep going with this KO only if there is at least one replicate
  if (replicatesPerKO > 1) {
  
  # add gene id as a column to merge tables
  ko.geneExp$gene_id = rownames(ko.geneExp)
  
  #join datasets
  geneExp = left_join(geneNmes.df, control.geneExp, by = "gene_id") %>%
    left_join(., ko.geneExp, by = "gene_id")
  
  #fid the column names
  colnames(geneExp)[1] = "NAME"
  
  write.table(geneExp, file = paste("03_extractedData/190321_KO-RNAseq/GSEA/geneExpressionFiles/", KO, "_gct.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")  
  
  ####now make the cls file####
  
  #detemrine number of samples
  numSamplesInFile = 10 + replicatesPerKO 
  
  #Create empty file
  clsFile = matrix(NA, 3, numSamplesInFile) #create file
  clsFile = data.frame(clsFile) #convert to data frame
  clsFile[,] = clsFile[is.na(clsFile[,])] = "" #make the data frame empty
  
  #fill in values into data frame
  clsFile[1,1] = numSamplesInFile 
  clsFile[1,2] = 2
  clsFile[1,3] = 1
  
  clsFile[2,1] = "#"
  clsFile[2,2] = "wt"
  clsFile[2,3] = "ko"
  
  clsFile[3,1:10] = "wt"
  clsFile[3,11:numSamplesInFile] = "ko"
  
  #save file
  write.table(clsFile, file = paste("03_extractedData/190321_KO-RNAseq/GSEA/phenotypeFiles/",KO, "_phenotype.cls", sep = ""), sep = " ", quote = F, col.names = F, row.names = F)

  }
  
}


