##The goal of this script is to perform GSEA in each of the KOs using the NGFR clusters as the gene sets

#Load groupings
groupings = read.table("02_metadata/190408_RNAseq_otherSamples/RNAseqGroupings.txt", header = T)
groupings$sample = as.character(groupings$sample)
groupings$baseline = as.character(groupings$baseline)

#Load metadata table
metadata = read.table("02_metadata/190408_RNAseq_otherSamples/shaffer_samplesSelected.txt", header = T)
metadata$sampleID = as.character(metadata$sampleID)
metadata$sampleID_new = as.character(metadata$sampleID_new)

#Load gene expresion data
geneExpData = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/normalizedCounts/metedData_normalized.txt", header = T)
geneExpData$gene_ID = as.character(geneExpData$gene_ID)
geneExpData$sampleID = as.character(geneExpData$sampleID)
geneExpData$counts_rpm = as.numeric(as.character(geneExpData$counts_rpm))

#Add new sample IDs to gene expression data
geneExpData = inner_join(geneExpData, metadata, by = "sampleID")

#Load name conversion file
nameConversionKey = read.table("02_metadata/190408_RNAseq_otherSamples/hg19_nameConvertionFile.txt", header = T)
nameConversionKey$gene_id = as.character(nameConversionKey$gene_id)

#Add new names to gene expression data
geneExpData = left_join(nameConversionKey, geneExpData, by = c("gene_id" = "gene_ID"))

#spread the table
geneExpData = geneExpData %>%
  dplyr::select(gene_id, GeneSymbol, sampleID_new, counts_rpm) %>%
  distinct()
geneExpData = geneExpData %>%
  spread(sampleID_new, counts_rpm)

#Clean up the table
geneExpData[is.na(geneExpData)] = 0 #convert NAs to 0
geneExpData = geneExpData[rowSums(geneExpData[4:ncol(geneExpData)]) > 0, ] #keep only genes where the sum of counts across all samples is at least #

# This tells me which genes remain duplicated.
tmp = geneExpData %>%
  dplyr::select(GeneSymbol) %>%
  mutate(dup = duplicated(GeneSymbol)) %>%
  filter(dup == T)
removeDuplicatedGene = tmp$GeneSymbol

#Remove the duplicated gene
geneExpData = geneExpData %>%
  filter(!GeneSymbol %in% removeDuplicatedGene) %>%
  distinct()

#Now convert the gene Symbol into the row names and eliminate unnecessary fields
row.names(geneExpData) = geneExpData$GeneSymbol #add rownames to prevent losing them
geneExpData$gene_id = NULL #eliminate column
geneExpData$GeneSymbol = NULL #eliminate column


#Make list of samples to loop through
RNAseqSamples = groupings$sample

for(currentSample in RNAseqSamples) {
  
  # #for testing
  # currentSample = RNAseqSamples[1]
  
  #create the first two columns of the output file
  geneNmes.df = data.frame(rownames(geneExpData))
  colnames(geneNmes.df) = "gene_id"
  
  geneNmes.df = geneNmes.df %>%
    mutate(DESCRIPTION = "na")
  
  controlGroup = groupings %>%
    filter(sample == currentSample)
  
  #create the table of controls
  control.geneExp = geneExpData %>%
    dplyr::select(contains(controlGroup$baseline))
  
  control.geneExp = control.geneExp[, colSums(control.geneExp) > 0]
  replicatesPerControls = ncol(control.geneExp)
  control.geneExp$gene_id = rownames(control.geneExp)
  
  #Select columns of interest
  ko.geneExp = geneExpData %>%
    dplyr::select(contains(controlGroup$sample))
  
  #Keep only columns that have content
  ko.geneExp = ko.geneExp[, colSums(ko.geneExp) > 0]
  ko.geneExp = data.frame(ko.geneExp)
  
  #Determines the number of biological replicates in a given KO
  replicatesPerKO = ncol(ko.geneExp)  #Im eliminating samples where I dont hve a replicate. 
  
  ####
  
  #Make GCT file and save it. BUT, keep going with this KO only if there is at least one replicate
  if (replicatesPerKO > 1) {
    
    # add gene id as a column to merge tables
    ko.geneExp$gene_id = rownames(ko.geneExp)
    
    #join datasets
    geneExp = left_join(geneNmes.df, control.geneExp, by = "gene_id") %>%
      left_join(., ko.geneExp, by = "gene_id")
    
    #fid the column names
    colnames(geneExp)[1] = "NAME"
    
    write.table(geneExp, file = paste("03_extractedData/190408_RNAseq_otherSamples/otherDatatsets/GSEA/geneExpressionFiles/", currentSample ,"_gct.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")  
    
    ####now make the cls file####
    
    #detemrine number of samples
    numSamplesInFile = replicatesPerControls + replicatesPerKO 
    
    #Create empty file
    clsFile = matrix(NA, 3, numSamplesInFile) #create file
    clsFile = data.frame(clsFile) #convert to data frame
    clsFile[,] = clsFile[is.na(clsFile[,])] = "" #make the data frame empty
    
    #fill in values into data frame
    clsFile[1,1] = numSamplesInFile 
    clsFile[1,2] = 2
    clsFile[1,3] = 1
    
    clsFile[2,1] = "#"
    clsFile[2,2] = "control"
    clsFile[2,3] = "condition"
    
    clsFile[3,1:replicatesPerControls] = "control"
    clsFile[3,(replicatesPerControls + 1):numSamplesInFile] = "condition"
    
    #save file
    write.table(clsFile, file = paste("03_extractedData/190408_RNAseq_otherSamples/otherDatatsets/GSEA/phenotypeFiles/", currentSample ,"_phenotype.cls", sep = ""), sep = " ", quote = F, col.names = F, row.names = F)
    
  }
  
  
}





