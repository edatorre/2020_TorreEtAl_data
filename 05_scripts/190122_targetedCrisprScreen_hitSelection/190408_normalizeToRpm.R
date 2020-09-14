#The goal of this script if to normalize data to reads per million

#Load data
rawReads = read.table("01_rawData/190122_targetedCrisprScreen_hitSelection/rawCounts_modifiedNames/rawReads_combined.txt", header = T)

#Make tall table
rawReads_tall = rawReads %>%
  gather("sample", "counts", 2:ncol(rawReads))

#Obtain total number of reads per sample
rawReads_tall = rawReads_tall %>%
  group_by(sample) %>%
  mutate(totalReads = sum(counts)) %>%
  ungroup(sample)

#Normalize to rpm
rawReads_tall = rawReads_tall %>%
  mutate(counts_rpm = counts * (10^6 / totalReads))

#convert back to matrix form
rawReads_tall = rawReads_tall %>%
  dplyr::select(-counts, -totalReads) %>%
  spread(sample, counts_rpm)

#save table
write.table(rawReads_tall, "01_rawData/190122_targetedCrisprScreen_hitSelection/rawCounts_modifiedNames/normalizedReads_combined.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
