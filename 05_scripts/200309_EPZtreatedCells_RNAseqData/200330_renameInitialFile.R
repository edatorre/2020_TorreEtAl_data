#The goal of this script is to convert sample names

#Read in raw counts
allCounts = read.table(file = "01_rawData/200309_EPZtreatedCells_RNAseqData/counts/meltedData.tsv", header = TRUE, sep = "\t")

#Read in the sample name conversion sheet
nameKey = read.table(file = "02_metadata/200309_EPZtreatedCells_RNAseqData/200312_listOfSamples.txt", header = TRUE, sep = "\t")

#Join Tables
allCounts = left_join(allCounts, nameKey, by = c("sampleID" = "sample" ))

#Remove unnecessary columns
allCounts = allCounts %>%
  dplyr::select(-experiment, -sampleID, -sampleType)

#rename column
colnames(allCounts)[3] = "sampleID"

#save table
write.table(allCounts, file = "01_rawData/200309_EPZtreatedCells_RNAseqData/counts/meltedData_renamed.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
