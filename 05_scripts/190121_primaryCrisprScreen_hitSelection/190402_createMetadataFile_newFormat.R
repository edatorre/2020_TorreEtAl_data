# #Create metadata file
# 
# rm(list = ls())
# 
# library(tidyverse)

#Define the positive controls
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")

#Load screen results for all tagrets: NGFRhi EGFRhi screen
stateScreen_epi = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_epigenetic_allInfo.txt", header = T, sep = "\t")
stateScreen_kinase = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_kinase_allInfo.txt", header = T, sep = "\t")
stateScreen_tf = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_TF_allInfo.txt", header = T, sep = "\t")
stateScreen = bind_rows(stateScreen_epi, stateScreen_kinase, stateScreen_tf)

stateScreen = stateScreen %>%
  dplyr::select(domain, library, medianlFC) %>%
  distinct()
colnames(stateScreen) = c("domain", "library","medianlFC_DP")

stateScreen = stateScreen %>%
  filter(!domain %in% positiveControls) %>%
  filter(!domain == "neg")

#Load screen results for all tagrets: Vemurafenib screen
fateScreen_epi = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_epigenetic_allInfo.txt", header = T, sep = "\t")
fateScreen_kinase = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_kinase_allInfo.txt", header = T, sep = "\t")
fateScreen_tf = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_TF_allInfo.txt", header = T, sep = "\t")
fateScreen = bind_rows(fateScreen_epi, fateScreen_kinase, fateScreen_tf)

fateScreen = fateScreen %>%
  dplyr::select(domain, library, medianlFC) %>%
  distinct()
colnames(fateScreen) = c("domain", "library","medianlFC_VemR")

fateScreen = fateScreen %>%
  filter(!domain %in% positiveControls) %>%
  filter(!domain == "neg")

screenResults = left_join(stateScreen, fateScreen, by = c("domain", "library"))

#Load sliding scale tables: DP
DP_epigenetic_table = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DP_epigeneticLibrary_slidingScaleTable.txt", header = T)
DP_kinase_table = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DP_kinaseLibrary_slidingScaleTable.txt", header = T)
DP_tf_table = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DP_TFLibrary_slidingScaleTable.txt", header = T)

#combine tables
DP_slidingScale = bind_rows(DP_epigenetic_table, DP_kinase_table, DP_tf_table)
DP_slidingScale = DP_slidingScale %>%
  dplyr::select(-gRNA) %>%
  distinct()
colnames(DP_slidingScale) = c("domain", "library", "DP_hitStatus_050", "DP_hitStatus_066", "DP_hitStatus_075")

#Load sliding scale tables: VemR
VemR_epigenetic_table = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemR_epigeneticLibrary_slidingScaleTable.txt", header = T)
VemR_kinase_table = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemR_kinaseLibrary_slidingScaleTable.txt", header = T)
VemR_tf_table = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemR_TFLibrary_slidingScaleTable.txt", header = T)

#combine tables
VemR_slidingScale = bind_rows(VemR_epigenetic_table, VemR_kinase_table, VemR_tf_table)
VemR_slidingScale = VemR_slidingScale %>%
  dplyr::select(-gRNA) %>%
  distinct()
colnames(VemR_slidingScale) = c("domain", "library", "VemR_hitStatus_050", "VemR_hitStatus_066", "VemR_hitStatus_075")

#combine tables
screenResults = left_join(screenResults, DP_slidingScale, by = c("domain", "library"))
screenResults = left_join(screenResults, VemR_slidingScale, by = c("domain", "library"))

#Add directionality info to hits
screenResults = screenResults %>%
  mutate(EffectOnScreen_DP = ifelse(medianlFC_DP < 0, "down", "up")) %>%
  mutate(EffectOnScreen_VemR = ifelse(medianlFC_VemR < 0, "down", "up"))

#Add tier information
DP_tier75 = screenResults %>%
  filter(DP_hitStatus_075 == "yes") %>%
  dplyr::select(domain) %>%
  mutate(DPtier = "tier075")

DP_tier66 = screenResults %>%
  filter(DP_hitStatus_075 == "No" & DP_hitStatus_066 == "yes") %>%
  dplyr::select(domain) %>%
  mutate(DPtier = "tier066")

DP_tier50 = screenResults %>%
  filter(DP_hitStatus_066 == "No" & DP_hitStatus_050 == "yes") %>%
  dplyr::select(domain) %>%
  mutate(DPtier = "tier050")

DP_tierNotHit = screenResults %>%
  filter(DP_hitStatus_050 == "No") %>%
  dplyr::select(domain) %>%
  mutate(DPtier = "tierNotHit")

DP_tiers = bind_rows(DP_tier75, DP_tier66, DP_tier50, DP_tierNotHit)

VemR_tier75 = screenResults %>%
  filter(VemR_hitStatus_075 == "yes") %>%
  dplyr::select(domain) %>%
  mutate(VemRtier = "tier075")

VemR_tier66 = screenResults %>%
  filter(VemR_hitStatus_075 == "No" & VemR_hitStatus_066 == "yes") %>%
  dplyr::select(domain) %>%
  mutate(VemRtier = "tier066")

VemR_tier50 = screenResults %>%
  filter(VemR_hitStatus_066 == "No" & VemR_hitStatus_050 == "yes") %>%
  dplyr::select(domain) %>%
  mutate(VemRtier = "tier050")

VemR_tierNotHit = screenResults %>%
  filter(VemR_hitStatus_050 == "No") %>%
  dplyr::select(domain) %>%
  mutate(VemRtier = "tierNotHit")

VemR_tiers = bind_rows(VemR_tier75, VemR_tier66, VemR_tier50, VemR_tierNotHit)

#Add tier info to main file
screenResults = left_join(screenResults, DP_tiers, by = "domain")
screenResults = left_join(screenResults, VemR_tiers, by = "domain")

#Fix a few of the names
screenResults$domain = gsub("BRD2_BD1", "BRD2-BD1", screenResults$domain)
screenResults$domain = gsub("BRD2_BD2", "BRD2-BD2", screenResults$domain)
screenResults$domain = sapply(strsplit(screenResults$domain,"_"), `[`, 1)
colnames(screenResults)[1] = "geneName"

#save file as metadata file
write.table(screenResults, "02_metadata/190402_metadataOfHits.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

