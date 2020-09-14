# #The goal of this script is to analyse NGFR IF data.

#Define the positive and negative controls
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
negativeControls = c("WT", "neg")

#Load name conversion file -> this file helps interconvert different version of the gRNA names
nameConverstion = read.table("02_metadata/190402_screenTargets_geneDomainTarget_nameConversion.txt", header = T)

#Load raw Immunofluorescence data
plate1 = read.csv("01_rawData/180925_NGFR-IF/plate1/myTable.csv")
plate2 = read.csv("01_rawData/180925_NGFR-IF/plate2/myTable.csv")
plate3 = read.csv("01_rawData/180925_NGFR-IF/plate3/myTable.csv")

#add plate ID
plate1 = plate1 %>%
  mutate(plateID = "plate1")
plate2 = plate2 %>%
  mutate(plateID = "plate2")
plate3 = plate3 %>%
  mutate(plateID = "plate3")

#combine data
allData = bind_rows(plate1, plate2, plate3)

#Load metadata
plateMaps = read.table("02_metadata/180925_NGFR-IF/NGFR-IF_plateMaps.txt", header = T)
nameKey = read.table("02_metadata/180925_NGFR-IF/NGFR-IF_sampleKey.txt", header = T)
metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)

#Fix names and add name conversions
colnames(nameKey)[3] = "gene"
nameKey = left_join(nameKey, nameConverstion, by = c("gRNA", "gene"))

plateMaps = left_join(plateMaps, nameKey, by = "sampleID")
allData = left_join(plateMaps, allData, by = c("plateID", "WellNumber" = "wellNumber"))

#add number of cells per well
allData = allData %>%
  group_by(gRNA) %>%
  mutate(numCells = n()) %>%
  #ungroup(sampleName) %>%
  ungroup(gRNA)

#Obtain summary of cells per well. 
cellsPerGuide_summary = allData %>%
  dplyr::select(target, gRNA, numCells) %>%
  distinct()

#Filter out samples with low number numbe rof cells
allData = allData %>%
  filter(numCells >= 500)

#Separate negative controls & detemrine number of cells that represents 1%
negControls = allData %>%
  filter(target == "neg") %>%
  group_by(gRNA) %>%
  mutate(top1 = round((numCells - (0.01 * numCells)))) %>%
  mutate(rankings = rank(ngfr_IFquant_intensities)) %>%
  filter(rankings >= top1) %>%
  ungroup(gRNA)

#Select the lowest immunofluorescence intensity value in the filtered ranked samples (the top1%)
selectMinValue = negControls %>%
  group_by(gRNA) %>%
  summarize(min(ngfr_IFquant_intensities)) %>%
  ungroup(gRNA)

#Plot histogram of min NGFR intensities to determine if there are outliers whose minimum intesity value is very different from the rest of the populaiton of controls
controlIntensities = ggplot(data = selectMinValue, aes(x = selectMinValue$`min(ngfr_IFquant_intensities)`)) +
  geom_histogram(bins = 10) +
  xlim(0,2500) + 
  theme_classic()
controlIntensities

#obtain a median across controls
selectMinValue = selectMinValue %>%
  mutate(IFthreshold = median(`min(ngfr_IFquant_intensities)`)) %>%
  dplyr::select(IFthreshold) %>%
  distinct()

#Keep only cells above the theshold.
cellAboveThreshold = allData %>%
  filter(ngfr_IFquant_intensities >= selectMinValue$IFthreshold)

#Calculate # of NGFR-hi cells per sample
cellAboveThreshold = cellAboveThreshold %>%
  group_by(gRNA) %>%
  mutate(numNGFRhiCells = n()) %>%
  ungroup(gRNA)

#Convert # NGFR-hi cells to a frequency
cellAboveThreshold = cellAboveThreshold %>%
  mutate(freqNGFRhi = numNGFRhiCells / numCells) %>%
  dplyr::select(gRNA, target, numCells, numNGFRhiCells, freqNGFRhi) %>%
  distinct()

#Obtain mean minimum intensity across negative controls
FractionNGFRhhiOnControls = cellAboveThreshold %>%
  filter(target == "neg") %>%
  summarise(mean(freqNGFRhi))

#Obtain fold change for all samples
cellAboveThreshold = cellAboveThreshold %>%
  mutate(log2_NGFRhi_FC = log2(freqNGFRhi / FractionNGFRhhiOnControls$`mean(freqNGFRhi)`))

#Obtain median and SD for all targets
cellAboveThreshold = cellAboveThreshold %>%
  group_by(target) %>%
  mutate(nReplicates = n_distinct(gRNA)) %>%
  mutate(meanNumCells = mean(numCells)) %>%
  mutate(meanlFC = mean(log2_NGFRhi_FC)) %>%
  mutate(sdlFC = sd(log2_NGFRhi_FC)) %>%
  mutate(selFC = sdlFC / sqrt(nReplicates)) %>%
  ungroup(target)

cellAboveThreshold = left_join(cellAboveThreshold, metadata, by = c("target" = "geneName"))
cellAboveThreshold$gRNA = as.character(cellAboveThreshold$gRNA)
cellAboveThreshold$target = as.character(cellAboveThreshold$target)

#sumarizedResults
sumarizedResults = cellAboveThreshold %>%
  dplyr::select(target, gRNA, log2_NGFRhi_FC, meanNumCells, meanlFC, sdlFC, selFC) %>% 
  distinct()

#remove positive controls
sumarizedResults = sumarizedResults %>%
  filter(!target %in% positiveControls) 
  #filter(!target %in% negativeControls)

#sumarizedResults
sumarizedResults_targetsOnly = cellAboveThreshold %>%
  dplyr::select(target, meanNumCells, meanlFC, selFC) %>%
  distinct()

#remove positive controls
sumarizedResults_targetsOnly = sumarizedResults %>%
  dplyr::select(target, meanlFC, sdlFC, selFC) %>%
  distinct()

#Save table
write.table(cellAboveThreshold, file = "03_extractedData/180925_NGFR-IF/allResults.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(sumarizedResults, file = "03_extractedData/180925_NGFR-IF/sumarizedResults.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(sumarizedResults_targetsOnly, file = "03_extractedData/180925_NGFR-IF/sumarizedResults_targetsOnly.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#merge sumarized results with metadata
sumarizedResults = left_join(sumarizedResults, metadata, by = c("target" = "geneName"))

####Plot results####
#Clean up a table 
sumarizedResults$sdlFC[is.na(sumarizedResults$sdlFC)] = 0
sumarizedResults = sumarizedResults %>%
  dplyr::select(target, gRNA, log2_NGFRhi_FC, meanlFC, sdlFC, selFC, EffectOnScreen_DP, EffectOnScreen_VemR, DP_hitStatus_075, DP_hitStatus_066, DP_hitStatus_050, VemR_hitStatus_075, VemR_hitStatus_066, VemR_hitStatus_050)

#PBRM1 and EP300 repeat several times in the table because one of the domains shows one effect and another ones shows a different effect. This is just an issue for organization but does not affect the IF results. I will remove one of avoid data repetition on the graph.
sumarizedResults = sumarizedResults %>%
  filter(!(target == "PBRM1" & EffectOnScreen_VemR == "down")) %>%
  filter(!(target == "EP300" & DP_hitStatus_075 == "yes"))
sumarizedResults = sumarizedResults %>%
  distinct()

#### Make plot organized by DP screen results ####
#make plotting table for the selected screen: DP
plottingTable_DP_tier1 = sumarizedResults %>%
  filter(DP_hitStatus_075 == "yes") %>%
  mutate(tier_DP = "tier1")
break1 = nrow(plottingTable_DP_tier1)

plottingTable_DP_tier2 = sumarizedResults %>%
  filter(DP_hitStatus_075 == "No" & DP_hitStatus_066 == "yes") %>%
  filter(!target %in% plottingTable_DP_tier1$target) %>%
  mutate(tier_DP = "tier2")
break2 = nrow(plottingTable_DP_tier2)

plottingTable_DP_tier3 = sumarizedResults %>%
  filter(DP_hitStatus_066 == "No" & DP_hitStatus_050 == "yes") %>%
  filter(!target %in% plottingTable_DP_tier2$target) %>%
  mutate(tier_DP = "tier3")
break3 = nrow(plottingTable_DP_tier3)

plottingTable_DP_tier4 = sumarizedResults %>%
  filter(DP_hitStatus_050 == "No") %>%
  filter(!target %in% plottingTable_DP_tier3$target) %>%
  mutate(tier_DP = "tier4")
break4 = nrow(plottingTable_DP_tier4)

plottingTable_DP_tierControl = sumarizedResults %>%
  filter(target == "neg" | target == "WT") %>%
  mutate(tier_DP = "xControl")

#Combine rows with new tier labels
plottingTable_DP = bind_rows(plottingTable_DP_tier1, plottingTable_DP_tier2, plottingTable_DP_tier3, plottingTable_DP_tier4, plottingTable_DP_tierControl)
plottingTable_DP = plottingTable_DP %>%
  distinct()

#Create tmp table to organize targets for the plot
tmp = plottingTable_DP %>%
  dplyr::select(-gRNA, -log2_NGFRhi_FC) %>%
  distinct() %>%
  group_by(tier_DP) %>%
  arrange(EffectOnScreen_DP, target) %>%
  mutate(plottingID = 1:(n_distinct(target))) %>%
  ungroup()

tmp = tmp %>%
  dplyr::select(target, plottingID)

plottingTable_DP = left_join(plottingTable_DP, tmp, by = "target")


#Make plot with all data points
NGFRbarplot = ggplot() +
  geom_point(data = plottingTable_DP, aes(x = plottingID, y = log2_NGFRhi_FC, color = EffectOnScreen_DP), stat = "identity", size = 0.5) +
  geom_text(data = plottingTable_DP, aes(x = plottingID, y = meanlFC, label = target, color = EffectOnScreen_DP), angle = 90, size = 2) +
  facet_grid(. ~ tier_DP) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_classic()
NGFRbarplot
ggsave("04_plots/180925_NGFR-IF/190403_NGFR-IF_allKOs_byDP_byTiers_onlyPoints.PDF", height = 3, width = 14)

#Make plot with bars and superimposed individual data points
tmp2 = plottingTable_DP %>%
  dplyr::select(-gRNA, -log2_NGFRhi_FC) %>%
  distinct()

NGFRbarplot = ggplot(data = tmp2, aes(x = plottingID, y = meanlFC, fill = EffectOnScreen_DP)) +
  geom_bar(stat = "identity") +
  geom_text(data = tmp2, aes(x = plottingID, y = meanlFC, label = target), angle = 90) +
  geom_point(data = plottingTable_DP, aes(x = plottingID, y = log2_NGFRhi_FC), stat = "identity", size = 0.5) +
  facet_grid(. ~ tier_DP) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_classic()
NGFRbarplot
ggsave("04_plots/180925_NGFR-IF/190403_NGFR-IF_allKOs_byDP_byTiers_pointsAndBars.PDF", height = 3, width = 14)

#### Make plot organized by VemR screen results ####
#make plotting table for the selected screen: VemR
plottingTable_VemR_tier1 = sumarizedResults %>%
  filter(VemR_hitStatus_075 == "yes") %>%
  mutate(tier_VemR = "tier1")
break1 = nrow(plottingTable_VemR_tier1)

plottingTable_VemR_tier2 = sumarizedResults %>%
  filter(VemR_hitStatus_075 == "No" & VemR_hitStatus_066 == "yes") %>%
  filter(!target %in% plottingTable_VemR_tier1$target) %>%
  mutate(tier_VemR = "tier2")
break2 = nrow(plottingTable_VemR_tier2)

plottingTable_VemR_tier3 = sumarizedResults %>%
  filter(VemR_hitStatus_066 == "No" & VemR_hitStatus_050 == "yes") %>%
  filter(!target %in% plottingTable_VemR_tier2$target) %>%
  mutate(tier_VemR = "tier3")
break3 = nrow(plottingTable_VemR_tier3)

plottingTable_VemR_tier4 = sumarizedResults %>%
  filter(VemR_hitStatus_050 == "No") %>%
  filter(!target %in% plottingTable_VemR_tier3$target) %>%
  mutate(tier_VemR = "tier4")
break4 = nrow(plottingTable_VemR_tier4)

plottingTable_VemR_tierControl = sumarizedResults %>%
  filter(target == "neg" | target == "WT") %>%
  mutate(tier_DP = "xControl")

plottingTable_VemR = bind_rows(plottingTable_VemR_tier1, plottingTable_VemR_tier2, plottingTable_VemR_tier3, plottingTable_VemR_tier4, plottingTable_VemR_tierControl)
plottingTable_VemR = plottingTable_VemR %>%
  distinct()


#Make temp table to organize hits for the plot
tmp3 = plottingTable_VemR %>%
  dplyr::select(-gRNA, -log2_NGFRhi_FC) %>%
  distinct() %>%
  group_by(tier_VemR) %>%
  arrange(EffectOnScreen_VemR, target) %>%
  mutate(plottingID = 1:(n_distinct(target))) %>%
  ungroup()

tmp3 = tmp3 %>%
  dplyr::select(target, plottingID)

plottingTable_VemR = left_join(plottingTable_VemR, tmp3, by = "target")

#Make plot
NGFRbarplot = ggplot() +
  geom_point(data = plottingTable_VemR, aes(x = plottingID, y = log2_NGFRhi_FC, color = EffectOnScreen_VemR), stat = "identity", size = 0.5) +
  geom_text(data = plottingTable_VemR, aes(x = plottingID, y = meanlFC, label = target, color = EffectOnScreen_VemR), angle = 90, size = 2) +
  facet_grid(. ~ tier_VemR) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_classic()
NGFRbarplot
ggsave("04_plots/180925_NGFR-IF/190403_NGFR-IF_allKOs_byVemR_byTiers_onlyPoints.PDF", height = 3, width = 14)

#This one is for bars and dots
tmp4 = plottingTable_VemR %>%
  dplyr::select(-gRNA, -log2_NGFRhi_FC) %>%
  distinct()

NGFRbarplot = ggplot(data = tmp4, aes(x = plottingID, y = meanlFC, fill = EffectOnScreen_VemR)) +
  geom_bar(stat = "identity") +
  geom_text(data = tmp4, aes(x = plottingID, y = meanlFC, label = target), angle = 90) +
  geom_point(data = plottingTable_VemR, aes(x = plottingID, y = log2_NGFRhi_FC), stat = "identity", size = 0.5) +
  facet_grid(. ~ tier_VemR) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme_classic()
NGFRbarplot
ggsave("04_plots/180925_NGFR-IF/190403_NGFR-IF_allKOs_byVemR_byTiers_pointsAndBars.PDF", height = 3, width = 14)
