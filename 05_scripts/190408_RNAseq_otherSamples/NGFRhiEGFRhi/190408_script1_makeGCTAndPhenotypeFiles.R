##The goal of this script is to perform GSEA in each of the KOs using the NGFR clusters as the gene sets

#Load gene expresion data
geneExpData = read.table("01_rawData/190408_RNAseq_otherSamples/NGFRandEGFRhi/normalizedCounts/normalizedCounts.txt", header = T)

#Load name conversion file
nameConversionKey = read.table("02_metadata/190321_KO-RNAseq/hg19_nameConvertionFile.txt", header = T)

#Add new names to gen expression data
geneExpData = left_join(nameConversionKey, geneExpData, by = "gene_id")

#Clean up the table
geneExpData[is.na(geneExpData)] = 0 #convert NAs to 0
geneExpData = geneExpData[rowSums(geneExpData[4:ncol(geneExpData)]) > 1, ] #keep only genes expressed in at least one of the samples

# This tells me which genes remain duplicated.
tmp = geneExpData %>%
  dplyr::select(GeneSymbol) %>%
  mutate(dup = duplicated(GeneSymbol)) %>%
  filter(dup == T) %>%
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

#create the first two columns of the output file
geneNmes.df = data.frame(rownames(geneExpData))
colnames(geneNmes.df) = "gene_id"

geneNmes.df = geneNmes.df %>%
  mutate(DESCRIPTION = "na")

#create the table of controls
control.geneExp = geneExpData %>%
  dplyr::select(contains("neg"))

control.geneExp = control.geneExp[, colSums(control.geneExp) > 0]
control.geneExp$gene_id = rownames(control.geneExp)

  
#Select columns of interest
ko.geneExp = geneExpData %>%
  dplyr::select(contains("DP")) #DP stands for double positive (NGFRhi and EGFRhi)

#Keep only columns that have content
ko.geneExp = ko.geneExp[, colSums(ko.geneExp) > 0]
ko.geneExp = data.frame(ko.geneExp)

#Determines the number of biological replicates in a given KO
replicatesPerKO = ncol(ko.geneExp)  #Im eliminating samples where I dont have a replicate. 

#Make GCT file and save it. BUT, keep going with this KO only if there is at least one replicate
if (replicatesPerKO > 1) {

  # add gene id as a column to merge tables
  ko.geneExp$gene_id = rownames(ko.geneExp)
  
  #join datasets
  geneExp = left_join(geneNmes.df, control.geneExp, by = "gene_id") %>%
    left_join(., ko.geneExp, by = "gene_id")
  
  #fid the column names
  colnames(geneExp)[1] = "NAME"
  
  write.table(geneExp, file = paste("03_extractedData/190408_RNAseq_otherSamples/NGFRandEGFRhi/GSEA/geneExpressionFiles/NGFRhiEGFRhi_gct.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")  
  
  ####now make the cls file####
  
  #detemrine number of samples
  numSamplesInFile = 3 + replicatesPerKO 
  
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
  clsFile[2,3] = "NGFRhiEGFRhi"
  
  clsFile[3,1:3] = "wt"
  clsFile[3,4:numSamplesInFile] = "NGFRhiEGFRhi"
  
  #save file
  write.table(clsFile, file = paste("03_extractedData/190408_RNAseq_otherSamples/NGFRandEGFRhi/GSEA/phenotypeFiles/NGFRhiEGFRhi_phenotype.cls", sep = ""), sep = " ", quote = F, col.names = F, row.names = F)

}



