#Load plate map
plateMap = as.tibble(read.table("02_metadata/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/colonyGrowthMetadata.txt", header = T))
plateMap = plateMap %>%
  filter(Plate == "plate3")

#Read in baseline plate 
baselineData = as_tibble(read.csv("01_rawData/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/summarizedResults/baseline_well13-18_plateSummary.csv"))
baselineData = baselineData[,1:2]
colnames(baselineData) = c("well", "cellNumber_baseline")

#Read in colony growth data
colonyGrowthData = as_tibble(read.csv("01_rawData/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/summarizedResults/onDrug_well13-18_plateSummary.csv"))
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
write.table(mergedData, file ="03_extractedData/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/plate3_extractedData.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
