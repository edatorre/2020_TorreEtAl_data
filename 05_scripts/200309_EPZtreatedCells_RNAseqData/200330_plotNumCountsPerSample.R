# #The goal of this script is to plot the number of reads per library. 
# 
# rm(list=ls())
# 
# #Load libraries that will be used in the scripts.
# library(tidyverse)
# #set the working directory. 
# setwd("/Users/etorre/Dropbox (RajLab)/Eduardo_shared/et_melanoma_PlasticityAndReprograming/190227_CRISPR-KO_run1to6/")

#Read in raw counts and metadata
allCounts = read.table(file = "01_rawData/200309_EPZtreatedCells_RNAseqData/counts/meltedData_renamed.tsv", header = TRUE, sep = "\t")


#determine number of reads per sample
allCounts = allCounts %>%
  group_by(sampleID) %>%
  mutate(totalNumReads = sum(counts)) %>%
  dplyr::select(sampleID, totalNumReads) %>%
  distinct()

#Determine total number of counts
x_sum = sum(allCounts$totalNumReads, na.rm = TRUE)
x_mean = mean(allCounts$totalNumReads, na.rm = TRUE)

#Plot counts per library
plot01 = ggplot(data = allCounts, aes(x = sampleID, y = totalNumReads))+
  geom_bar(stat = "identity") +
  theme_classic()
plot(plot01)
ggsave(plot01, file = "04_plots/200309_EPZtreatedCells_RNAseqData/numReadsPerSample.PDF", height = 7, width = 7)

#save table of total number of reads per sample
write.table(allCounts, file = "01_rawData/200309_EPZtreatedCells_RNAseqData/counts/countsPerSample.txt", sep = "\t", row.names = FALSE, quote = FALSE)



