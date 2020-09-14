#The goal of this script is to merge melted data from multiple runs into a single file. 

#Read in raw counts
run1Counts = read.table(file = "01_rawData/190321_KO-RNAseq/rawCounts/meltedData_run1.tsv", header = TRUE, sep = "\t")
run2Counts = read.table(file = "01_rawData/190321_KO-RNAseq/rawCounts/meltedData_run2.tsv", header = TRUE, sep = "\t")
run3Counts = read.table(file = "01_rawData/190321_KO-RNAseq/rawCounts/meltedData_run3.tsv", header = TRUE, sep = "\t")
run4Counts = read.table(file = "01_rawData/190321_KO-RNAseq/rawCounts/meltedData_run4.tsv", header = TRUE, sep = "\t")
run5Counts = read.table(file = "01_rawData/190321_KO-RNAseq/rawCounts/meltedData_run5.tsv", header = TRUE, sep = "\t")
run6Counts = read.table(file = "01_rawData/190321_KO-RNAseq/rawCounts/meltedData_run6.tsv", header = TRUE, sep = "\t")

#Read in gene name conversion file
conversionKey = read.table(file = "02_metadata/190321_KO-RNAseq/hg19_nameConvertionFile.txt", header = TRUE, sep = "\t")

#Merge datasets
allCounts = left_join(run1Counts, run2Counts, by = c("sampleID", "gene_id")) %>%
  left_join(. , run3Counts, by = c("sampleID", "gene_id")) %>%
  left_join(. , run4Counts, by = c("sampleID", "gene_id")) %>%
  left_join(. , run5Counts, by = c("sampleID", "gene_id")) %>%
  left_join(. , run6Counts, by = c("sampleID", "gene_id")) 
  
#Sum counts from all sequencing runs
allCounts = allCounts %>%
  mutate(counts = counts.x + counts.y + counts.x.x + counts.y.y + counts.x.x.x + counts.y.y.y) %>%
  select(sampleID, gene_id, counts)

#Remove rows that do not correspond to genes
allCounts = filter(allCounts, grepl("ENSG", gene_id))

#swap gene_id for geneName
allCounts = left_join(allCounts, conversionKey, by = "gene_id")
allCounts = allCounts %>%
  select(-gene_id) %>%
  rename(gene_id = GeneSymbol_gene_id)


#saveRawCounts table
write.table(allCounts, file = "01_rawData/190321_KO-RNAseq/rawCounts/rawCounts_Runs1to6.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)
