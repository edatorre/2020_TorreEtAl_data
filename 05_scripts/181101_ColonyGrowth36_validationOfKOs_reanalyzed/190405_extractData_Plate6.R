# #The goal of this script is to analyze and plot colony gorwth data
# rm(list = ls())
# 
# #load library
# library(tidyverse)

#Load plate map
plateMap = as.tibble(read.table("02_metadata/181101_ColonyGrowth36_validationOfKOs_reanalyzed/plateMaps/plate6.txt", header = T))

#Read in baseline plate 
baselineData = as_tibble(read.csv("01_rawData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/summarizedResults/row6_baseline_plateSummary.csv"))
baselineData = baselineData[,1:2]
colnames(baselineData) = c("well", "cellNumber_baseline")

#Read in colony growth data
colonyGrowthData = as_tibble(read.csv("01_rawData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/summarizedResults/row6_onPLX4032_plateSummary.csv"))
colnames(colonyGrowthData) = c("well", "cellNumber_onDrug", "colonyNumber_onDrug", "numCellsInsideColonies_onDrug", "avgCellsPerColony_onDrug", "cellOutsideColonies_onDrug")

#merge data
mergedData = left_join(plateMap, baselineData, by = "well")
mergedData = left_join(mergedData, colonyGrowthData, by = "well")

#compute metrics of resistance
mergedData = mergedData %>%
  mutate(totalRcells_norm = cellNumber_onDrug / cellNumber_baseline) %>%
  mutate(Rcolonies_norm = colonyNumber_onDrug * 10000 / cellNumber_baseline) %>%
  mutate(survivingCells_norm = cellOutsideColonies_onDrug / cellNumber_baseline)

#save output table
setwd("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/")
write.table(mergedData, file ="plate6_extractedData.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


