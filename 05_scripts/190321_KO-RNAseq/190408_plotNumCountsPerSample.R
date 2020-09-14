# #The goal of this script is to plot the number of reads per library. 


#Read in raw counts and metadata
allCounts = read.table(file = "01_rawData/190321_KO-RNAseq/rawCounts/rawCounts_Runs1to6.tsv", header = TRUE, sep = "\t")
metadata = read.table(file = "02_metadata/190402_metadataOfHits.txt", header = TRUE, sep = "\t")

#determine number of reads per sample
allCounts = allCounts %>%
  group_by(sampleID) %>%
  mutate(totalNumReads = sum(counts)) %>%
  select(sampleID, totalNumReads) %>%
  distinct()

#Determine total number of counts
x_sum = sum(allCounts$totalNumReads, na.rm = TRUE)
x_mean = mean(allCounts$totalNumReads, na.rm = TRUE)

#Plot counts per library
plot01 = ggplot(data = allCounts, aes(x = sampleID, y = totalNumReads))+
  geom_bar(stat = "identity") +
  theme_classic()
plot(plot01)
ggsave(plot01, file = "04_plots/190321_KO-RNAseq/numReadsPerSample.PDF", height = 7, width = 7)

#save table of total number of reads per sample
write.table(allCounts, file = "01_rawData/190321_KO-RNAseq/rawCounts/countsPerSample.txt", sep = "\t", row.names = FALSE, quote = FALSE)



