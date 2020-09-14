#The goal of this script is to obtain the mean fold change of sgRNAs in the primary screen including only sgRNAs used in the targeted screen.

#Load sgRNAs used in the targeted screen
targetedScreen_targets = read.table("02_metadata/TargetedScreen_targets.txt", header = T)

#Load and combine primary screen results
primaryScreen_DP_epi = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_epigenetic_allInfo.txt", header = T, sep = "\t")
primaryScreen_DP_kinase = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_kinase_allInfo.txt", header = T, sep = "\t")
primaryScreen_DP_tf = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_TF_allInfo.txt", header = T, sep = "\t")

primaryScreen_DP = bind_rows(primaryScreen_DP_epi, primaryScreen_DP_kinase, primaryScreen_DP_tf)

primaryScreen_VemR_epi = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_epigenetic_allInfo.txt", header = T, sep = "\t")
primaryScreen_VemR_kinase = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_kinase_allInfo.txt", header = T, sep = "\t")
primaryScreen_VemR_tf = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_TF_allInfo.txt", header = T, sep = "\t")

primaryScreen_VemR = bind_rows(primaryScreen_VemR_epi, primaryScreen_VemR_kinase, primaryScreen_VemR_tf)

####Work on DP results

#Keep only sgRNAs used in targeted screen
primaryScreen_DP = left_join(targetedScreen_targets, primaryScreen_DP, by = "gRNA")
primaryScreen_VemR = left_join(targetedScreen_targets, primaryScreen_VemR, by = "gRNA")

# tmp = anti_join(primaryScreen_VemR, primaryScreen_DP, by = "gRNA")

#Calculate median log2 FC by domain
primaryScreen_DP = primaryScreen_DP %>%
  group_by(domain) %>%
  mutate(medianlFC_forTargetedScreen = median(lFC)) %>%
  mutate(sdlFC_forTargetedScreen = sd(lFC)) %>%
  ungroup(domain)

primaryScreen_VemR = primaryScreen_VemR %>%
  group_by(domain) %>%
  mutate(medianlFC_forTargetedScreen = median(lFC)) %>%
  mutate(sdlFC_forTargetedScreen = sd(lFC)) %>%
  ungroup(domain)

#Limit table to columns to use later
primaryScreen_DP = primaryScreen_DP %>%
  dplyr::select(target, gene, domain, lFC, medianlFC, medianlFC_forTargetedScreen, sdlFC_forTargetedScreen)

primaryScreen_VemR = primaryScreen_VemR %>%
  dplyr::select(target, gene, domain, lFC, medianlFC, medianlFC_forTargetedScreen, sdlFC_forTargetedScreen)

#save tables
write.table(primaryScreen_DP, file = "03_extractedData/190121_primaryCrisprScreen_hitSelection/DPscreen_medianlFC_targetedScreenGuides.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(primaryScreen_VemR, file = "03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRscreen_medianlFC_targetedScreenGuides.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

