##The goal of this script is to normalize reads to rpm

#Load datasets
#This is the raw counts file
t01_allData = read.table("01_rawData/190408_RNAseq_otherSamples/NGFRandEGFRhi/rawCounts/meltedData.tsv", header = T)

#Count reads per library
t01_allData = t01_allData %>%
  group_by(sampleID) %>%
  mutate(totalReads = sum(counts)) %>%
  ungroup(sampleID)

#Normalize counts 
t01_allData = t01_allData %>%
  mutate(counts_rpm = (counts * 10^6) / totalReads )

#reformat file
t01_allData_matrix = t01_allData %>%
  dplyr::select(gene_id, sampleID,counts_rpm) %>%
  spread(sampleID, counts_rpm)

t01_allData_matrix = t01_allData_matrix[!grepl("_", t01_allData_matrix$gene_id),]

#Save table
write.table(t01_allData_matrix, file = "01_rawData/190408_RNAseq_otherSamples/NGFRandEGFRhi/normalizedCounts/normalizedCounts.txt", quote = FALSE, sep = "\t", col.names = T, row.names = FALSE)
