# #The goal of this script is to plot the colony growth data from colony growth 36.
# rm(list = ls())
# 
# library(tidyverse)

#Controls
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
negativeControls = c("WT", "neg", "Neg01", "Neg02")

#Load name conversion file -> this file helps interconvert different version of the gRNA names
nameConverstion = read.table("02_metadata/190402_screenTargets_geneDomainTarget_nameConversion.txt", header = T)

#load metadata
metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)  

#Add name conversion to metadata
metadata = left_join(metadata, nameConverstion, by = c("geneName" = "target"))

#load data
p1 = read.table("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/plate1_extractedData.txt", header = T)
p2 = read.table("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/plate2_extractedData.txt", header = T)
p3 = read.table("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/plate3_extractedData.txt", header = T)
p4 = read.table("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/plate4_extractedData.txt", header = T)
p5 = read.table("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/plate5_extractedData.txt", header = T)
p6 = read.table("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/plate6_extractedData.txt", header = T)

#Combine data
allData = bind_rows(p1, p2, p3, p4, p5, p6)

#Eliminate unnecessary field
allData = allData %>%
  dplyr::select(-well)

#obtain mean metrics of the negative controls
negativeControls = allData %>%
  filter(target == "Neg01" | target == "Neg02")

#summarized values of negative controls
meanTotalRcells_negControl = negativeControls %>%
  summarise(meanTotalRcells = mean(totalRcells_norm)) 

meanRcolonies_negControl = negativeControls %>%
  summarise(meanRcolonies = mean(Rcolonies_norm))

meanSurvivingCells_negControl = negativeControls %>%
  summarise(meanSurvivingCells = mean(survivingCells_norm))

#obtain log fold changes for all samples
allData = allData %>%
  mutate(totalRcells_lFC = log2(totalRcells_norm / meanTotalRcells_negControl$meanTotalRcells)) %>%
  mutate(Rcolonies_lFC = log2(Rcolonies_norm / meanRcolonies_negControl$meanRcolonies)) %>%
  mutate(survivingCells_lFC = log2(survivingCells_norm / meanSurvivingCells_negControl$meanSurvivingCells))

#Combine file with metadata
allData = left_join(allData, metadata, by = c("target" = "gene"))
allData = allData %>%
  dplyr::select(target, totalRcells_lFC, Rcolonies_lFC, survivingCells_lFC, EffectOnScreen_DP, EffectOnScreen_VemR, DPtier, VemRtier)
allData = allData %>%
  distinct()

##I am removing one of the copies of EP300 because its duplicated due to having multiple domains as targets that are hits. This does not affect the data or the plot.
allData = allData %>%
  filter(! (target == "EP300" & DPtier == "tier066"))

#Eliminate controls and samples not relevant in analysis (not hits form the screen). 
allData = allData %>%
  filter(!target %in% positiveControls)
  #filter(!target %in% negativeControls)

#Fix the names of the tiers to plot more easily
allData$DPtier = gsub("tier075", "tier1", allData$DPtier)
allData$DPtier = gsub("tier066", "tier2", allData$DPtier)
allData$DPtier = gsub("tier050", "tier3", allData$DPtier)
allData$DPtier = gsub("tierNotHit", "tier4", allData$DPtier)
allData$VemRtier = gsub("tier075", "tier1", allData$VemRtier)
allData$VemRtier = gsub("tier066", "tier2", allData$VemRtier)
allData$VemRtier = gsub("tier050", "tier3", allData$VemRtier)
allData$VemRtier = gsub("tierNotHit", "tier4", allData$VemRtier)

#save table in case you need to use it for something else
write.table(allData, file = "03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/colonyGrowthResults_allhits.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


#Add plotting IDs to the data
allData = allData %>%
  group_by(DPtier) %>%
  mutate(plottingID_DP = 1:n_distinct(target)) %>%
  ungroup(DPtier)

allData = allData %>%
  group_by(VemRtier) %>%
  mutate(plottingID_VemR = 1:n_distinct(target)) %>%
  ungroup(VemRtier)

#make plot: All hits by DP tiers
plotResistantColonies = ggplot() +
  geom_bar(data = allData, aes(x = plottingID_DP, y = Rcolonies_lFC, fill = EffectOnScreen_DP), stat = "identity") +
  geom_text(data = allData, aes(x = plottingID_DP, y = Rcolonies_lFC, label = target), size = 2, angle = 90) +
  facet_grid(. ~ DPtier) +
  # scale_fill_manual(values = c("control" = "gray", "down" = "black", "up" = "green", "." = "red")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = 0) + 
  theme_classic()
plotResistantColonies
ggsave("04_plots/181101_ColonyGrowth36_validationOfKOs_reanalyzed/colonyGrowthResults_DPhits_byTier_barplot.PDF", height = 2, width = 6)

#make plot: all hits by VemR tiers
plotResistantColonies = ggplot() +
  geom_bar(data = allData, aes(x = plottingID_VemR, y = Rcolonies_lFC, fill = EffectOnScreen_VemR), stat = "identity") +
  geom_text(data = allData, aes(x = plottingID_VemR, y = Rcolonies_lFC, label = target), size = 2, angle = 90) +
  facet_grid(. ~ VemRtier) +
  # scale_fill_manual(values = c("control" = "gray", "down" = "black", "up" = "green", "." = "red")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = 0) + 
  theme_classic()
plotResistantColonies
ggsave("04_plots/181101_ColonyGrowth36_validationOfKOs_reanalyzed/colonyGrowthResults_VemRhits_byTier_barplot.PDF", height = 2, width = 6)

