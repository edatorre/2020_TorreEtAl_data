# #The goal of this script is to merge melted data from multiple runs into a single file. 

#Read in raw counts
allCounts = read.table(file = "01_rawData/200309_EPZtreatedCells_RNAseqData/counts/meltedData_renamed.tsv", header = TRUE, sep = "\t")

#determine number of reads per sample
allCounts = allCounts %>%
  group_by(sampleID) %>%
  mutate(totalNumReads = sum(counts))

#normalize to rpm
allCounts = allCounts %>%
  mutate(counts_rpm = ((counts * 10^6)/totalNumReads))

#remove raw counts
allCounts_rpm = allCounts %>%
  dplyr::select(-counts, -totalNumReads)

#save rpm table
write.table(allCounts_rpm, file = "01_rawData/200309_EPZtreatedCells_RNAseqData/normalizedConts/rpmCounts.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#convert tall table to matrix
allCounts_rpmMatrix = allCounts_rpm %>%
  spread(sampleID, counts_rpm)

#save rpm matrix
write.table(allCounts_rpmMatrix, file = "01_rawData/200309_EPZtreatedCells_RNAseqData/normalizedConts/rpmCounts_matrix.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

