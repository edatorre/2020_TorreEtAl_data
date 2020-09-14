##The goal of this script is to normalize values to rpm


#Load raw counts
rawCounts = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/allCounts_meltedData.txt", header = T)

#Normalized to reads per million
rawCounts = rawCounts %>%
  group_by(sampleID) %>%
  mutate(totalReads = sum(counts)) %>%
  mutate(counts_rpm = (10^6 * counts) / totalReads) %>%
  ungroup(sampleID)

#remove unwanted column
rawCounts = rawCounts %>%
  dplyr::select(-counts, -totalReads)

rawCounts_matrix = rawCounts %>%
  spread(sampleID, counts_rpm)

#save file 
write.table(rawCounts, file = "01_rawData/190408_RNAseq_otherSamples/otherDatasets/normalizedCounts/metedData_normalized.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(rawCounts_matrix, file = "01_rawData/190408_RNAseq_otherSamples/otherDatasets/normalizedCounts/metedData_normalized_matrix.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
