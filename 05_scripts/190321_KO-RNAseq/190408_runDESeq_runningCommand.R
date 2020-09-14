#This script runs DEseq on all the samples

#Load metadata
t01_metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)

#Load list of targets
ListOfGenes = read.table("02_metadata/TargetedScreen_targets.txt", header = T)
ListOfGenes = ListOfGenes %>%
  dplyr::select(target) %>%
  distinct()

#Load raw counts
t01_allData = read.table("01_rawData/190321_KO-RNAseq/rawCounts/rawCounts_Runs1to6.tsv", header = T)

#Load counts per sample
countsPerSample = read.table("01_rawData/190321_KO-RNAseq/rawCounts/countsPerSample.txt", header = T)
countsPerSample$totalNumReads[is.na(countsPerSample$totalNumReads)] = 0

#Identify the samples with low number of reads
lowReadSamples = countsPerSample %>%
  filter(totalNumReads < 500000)
badSamples = as.vector(lowReadSamples$sampleID)

#Replace empty values with 0  
t01_allData[is.na(t01_allData)] <- 0

#Clean up and re-structure data.
#Convert tall table into matrix
allCounts_rawCountsMatrix = t01_allData %>%
  spread(sampleID, counts)

#Remove columns with bad samples: these samples had low reads
allCounts_rawCountsMatrix = allCounts_rawCountsMatrix %>%
  select(-one_of(badSamples))

#Make gene_id the rowname: prevents losing the gene_id during column selection
row.names(allCounts_rawCountsMatrix) = allCounts_rawCountsMatrix$gene_id
allCounts_rawCountsMatrix$gene_id = NULL

#Make Forloop to run DESeq
#Read in the function: runDESeq
source("05_scripts/190321_KO-RNAseq/190226_DEseq2.R")

#Run the function in each of the genes 
for (geneName in ListOfGenes$target) {
  
  print(paste("Running DEseq on ", geneName))
  flush.console()
  
  runDESeq(geneName)
  
}
