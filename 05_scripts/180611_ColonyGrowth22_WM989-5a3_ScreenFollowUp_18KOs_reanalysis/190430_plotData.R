
#Controls
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
negativeControls = c("WT", "neg", "Neg01", "Neg02", "EMX1")

#Load name conversion file -> this file helps interconvert different version of the gRNA names
nameConverstion = read.table("02_metadata/190402_screenTargets_geneDomainTarget_nameConversion.txt", header = T)

#load metadata
metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)

#Add name conversion to metadata
metadata = left_join(metadata, nameConverstion, by = c("geneName" = "target"))
metadata = metadata %>%
  dplyr::select(gene, EffectOnScreen_DP, EffectOnScreen_VemR, DPtier, VemRtier)
metadata = metadata %>%
  distinct()

#load data
p1 = read.table("03_extractedData/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/plate1_extractedData.txt", header = T)
p2 = read.table("03_extractedData/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/plate2_extractedData.txt", header = T)
p3 = read.table("03_extractedData/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/plate3_extractedData.txt", header = T)


#Combine data
allData = bind_rows(p1, p2, p3)

#Eliminate unnecessary field
allData = allData %>%
  dplyr::select(-well) %>%
  dplyr::select(-Well)

# #Fix a couple of the names to help them match the nomenclature in the dataset
# metadata$geneName = gsub("BRD2-BD1", "BRD2", metadata$geneName)

#obtain mean metrics of the negative controls
negativeControls = allData %>%
  filter(KO %in% negativeControls)

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
allData = left_join(allData, metadata, by = c("KO" = "gene"))

allData = allData %>%
  dplyr::select(KO, totalRcells_lFC, Rcolonies_lFC, survivingCells_lFC, EffectOnScreen_DP, EffectOnScreen_VemR, DPtier, VemRtier)
allData = allData %>%
  distinct()

#Eliminate controls and samples not relevant in analysis (not hits form the screen). 
allData = allData %>%
  filter(!KO %in% positiveControls)

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
write.table(allData, file = "03_extractedData/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/colonyGrowthResults_allhits.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


#Add plotting IDs to the data
allData = allData %>%
  group_by(DPtier) %>%
  arrange(EffectOnScreen_DP) %>%
  mutate(plottingID_DP = 1:n_distinct(KO)) %>%
  ungroup(DPtier)

allData = allData %>%
  group_by(VemRtier) %>%
  arrange(EffectOnScreen_VemR) %>%
  mutate(plottingID_VemR = 1:n_distinct(KO)) %>%
  ungroup(VemRtier)

#make plot: All hits by DP tiers
plotResistantColonies = ggplot() +
  geom_bar(data = allData, aes(x = plottingID_DP, y = Rcolonies_lFC, fill = EffectOnScreen_DP), stat = "identity") +
  geom_text(data = allData, aes(x = plottingID_DP, y = Rcolonies_lFC, label = KO), size = 2, angle = 90) +
  facet_grid(. ~ DPtier) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = 0) + 
  theme_classic()
plotResistantColonies
ggsave("04_plots/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/colonyGrowthResults_byDPTier_barplot_replicate1.PDF", height = 2, width = 6)

#make plot: all hits by VemR tiers
plotResistantColonies = ggplot() +
  geom_bar(data = allData, aes(x = plottingID_VemR, y = Rcolonies_lFC, fill = EffectOnScreen_VemR), stat = "identity") +
  geom_text(data = allData, aes(x = plottingID_VemR, y = Rcolonies_lFC, label = KO), size = 2, angle = 90) +
  facet_grid(. ~ VemRtier) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = 0) + 
  theme_classic()
plotResistantColonies
ggsave("04_plots/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/colonyGrowthResults_byVemRTier_barplot_replicate1.PDF", height = 2, width = 6)

