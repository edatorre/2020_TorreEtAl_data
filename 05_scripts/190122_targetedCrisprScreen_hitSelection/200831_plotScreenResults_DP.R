#The goal of this script is to plot the results from the targeted screen

#Load metadata
metadata = read.table("02_metadata/190402_metadataOfHits.txt", sep = "\t", header = T)

#Keep only metadata needed
metadata = metadata %>%
  filter(!(geneName == "EP300" & DPtier == "tier066")) %>%
  filter(!(geneName == "PBRM1" & DPtier == "tierNotHit")) %>%
  filter(!(geneName == "PBRM1" & library == "epigenetic"))

#Fix names of tiers for the plot
metadata$DPtier = gsub("tier075", "tier1_075", metadata$DPtier)
metadata$DPtier = gsub("tier066", "tier2_066", metadata$DPtier)
metadata$DPtier = gsub("tier050", "tier3_050", metadata$DPtier)
metadata$DPtier = gsub("tierNotHit", "tier4_NotHit", metadata$DPtier)

metadata$VemRtier = gsub("tier075", "tier1_075", metadata$VemRtier)
metadata$VemRtier = gsub("tier066", "tier2_066", metadata$VemRtier)
metadata$VemRtier = gsub("tier050", "tier3_050", metadata$VemRtier)
metadata$VemRtier = gsub("tierNotHit", "tier4_NotHit", metadata$VemRtier)

#Load results from primary screen and arrange table
primaryScreen = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPscreen_medianlFC_targetedScreenGuides.txt", header = T)
primaryScreen = primaryScreen %>%
  dplyr::select(target, lFC, medianlFC_forTargetedScreen, sdlFC_forTargetedScreen)
colnames(primaryScreen)[3] = "medianlFC"
colnames(primaryScreen)[4] = "sdlFC"
primaryScreen$cellLine = "primaryScreen"
primaryScreen = primaryScreen %>%
  distinct()

#Load results
cellLine1 = read.table("03_extractedData/190122_targetedCrisprScreen_hitSelection/DPanalysis_WM989.bulk_allInfo.txt", header = T)
cellLine1 = cellLine1 %>%
  dplyr::select(target, lFC, medianlFC, sdlFC) %>%
  distinct()
cellLine1$cellLine = "WM989.bulk"


#Load results
cellLine2 = read.table("03_extractedData/190122_targetedCrisprScreen_hitSelection/DPanalysis_X451Lu_allInfo.txt", header = T)
cellLine2 = cellLine2 %>%
  dplyr::select(target, lFC, medianlFC, sdlFC) %>%
  distinct()
cellLine2$cellLine = "X451Lu"

#Combine results and add metadata
screenResults = bind_rows(cellLine1, cellLine2, primaryScreen)
screenResults = left_join(screenResults, metadata, by = c("target" = "geneName"))

#Eliminate the lFC values to arrange the plot 
screenResults2 = screenResults %>%
  dplyr::select(-lFC) %>%
  distinct()

#Add plotting ID to all targets
tmp = screenResults2 %>%
  filter(cellLine == "X451Lu") %>%
  distinct() %>%
  group_by(DPtier) %>%
  arrange(EffectOnScreen_DP, target) %>%
  mutate(plottingID = 1:(n_distinct(target))) %>%
  ungroup(DPtier) %>%
  dplyr::select(target, plottingID)

screenResults2 = left_join(screenResults2, tmp, by = "target")
screenResults = left_join(screenResults, tmp, by = "target")


#PLot the results - bars and dots
plot_screenResults = ggplot(data = screenResults, aes(x = plottingID, y = medianlFC)) +
  geom_bar(data = screenResults, aes(x = plottingID, y = medianlFC, fill = cellLine), stat = "identity", position = position_dodge(1), alpha = 0.25) +
  geom_point(data = screenResults, aes(x = plottingID, y = lFC, fill = cellLine), stat = "identity", size = 0.25, position = position_dodge(1)) +
  geom_text(aes(x = plottingID, y = -3, label = target, color = EffectOnScreen_DP), size = 2) +
  scale_color_manual(values = c("down" = "#000000", "up" = "#00CC00", "NA" = "blue", "primaryScreen" = "#000000", "WM989.bulk" = "#f2a32c", "X451Lu" = "#2ce8f2")) +
  scale_fill_manual(values = c("primaryScreen" = "#000000", "WM989.bulk" = "#f2a32c", "X451Lu" = "#2ce8f2")) +
  geom_hline(yintercept = c(0), size = 0.25) +
  facet_grid(. ~ DPtier, scales = "free") +
  coord_flip() +
  theme_classic()
plot_screenResults
ggsave("04_plots/190122_targetedCrisprScreen_hitSelection/200831_DPscreen_twoCellLines.PDF", height = 8, width = 7)

#PLot the results - dots only
plot_screenResults = ggplot(data = screenResults, aes(x = plottingID, y = medianlFC)) +
  geom_point(data = screenResults, aes(x = plottingID, y = lFC, color = cellLine), stat = "identity", size = 0.25, position = position_dodge(1)) +
  geom_text(aes(x = plottingID, y = -3, label = target, color = EffectOnScreen_DP), size = 2) +
  scale_color_manual(values = c("down" = "#000000", "up" = "#00CC00", "NA" = "blue", "primaryScreen" = "#000000", "WM989.bulk" = "#f2a32c", "X451Lu" = "#2ce8f2")) +
  #scale_fill_manual(values = c("primaryScreen" = "#000000", "WM989.bulk" = "#f2a32c", "X451Lu" = "#2ce8f2")) +
  geom_hline(yintercept = c(0), size = 0.25) +
  facet_grid(. ~ DPtier, scales = "free") +
  coord_flip() +
  theme_classic()
plot_screenResults
ggsave("04_plots/190122_targetedCrisprScreen_hitSelection/200831_DPscreen_twoCellLines_dotsOnly.PDF", height = 8, width = 7)


