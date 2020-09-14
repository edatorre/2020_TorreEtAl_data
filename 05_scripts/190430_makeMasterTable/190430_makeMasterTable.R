#The goal of this script is to create a master table that contains most of the CRISPR results from the paper

#Load anme key
nameKey = read.table("02_metadata/screenTargets_geneDomainTarget_nameConversion.txt", header = T)

#load targets of primary screen
screenTargets = read.table("02_metadata/ScreenTargets_metadata.txt", header = T)

#Add name key to screen targets
masterTable = left_join(screenTargets, nameKey, by = c("gRNA", "gene", "domain"))

#Load primary screen results
primaryScreenResults_DP_epigeneticLibrary = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_epigenetic_allInfo.txt", header = T, sep = "\t")
primaryScreenResults_DP_kinaseLibrary = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_kinase_allInfo.txt", header = T, sep = "\t")
primaryScreenResults_DP_tfLibrary = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_TF_allInfo.txt", header = T, sep = "\t")

primaryScreenResults_DP = bind_rows(primaryScreenResults_DP_epigeneticLibrary, primaryScreenResults_DP_kinaseLibrary, primaryScreenResults_DP_tfLibrary)

primaryScreenResults_DP = primaryScreenResults_DP %>%
  dplyr::select(gRNA, library, lFC, medianlFC)
colnames(primaryScreenResults_DP) = c("gRNA", "library","primaryStateScreen_lFC", "primaryStateScreen_medianlFC")

primaryScreenResults_VemR_epigeneticLibrary = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_epigenetic_allInfo.txt", header = T, sep = "\t")
primaryScreenResults_VemR_kinaseLibrary = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_kinase_allInfo.txt", header = T, sep = "\t")
primaryScreenResults_VemR_tfLibrary = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_TF_allInfo.txt", header = T, sep = "\t")

primaryScreenResults_VemR = bind_rows(primaryScreenResults_VemR_epigeneticLibrary, primaryScreenResults_VemR_kinaseLibrary, primaryScreenResults_VemR_tfLibrary)

primaryScreenResults_VemR = primaryScreenResults_VemR %>%
  dplyr::select(gRNA, library, lFC, medianlFC)
colnames(primaryScreenResults_VemR) = c("gRNA", "library","primaryFateScreen_lFC", "primaryFateScreen_medianlFC")

#Add tiers to samples
primaryScreen_metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)
primaryScreen_metadata = primaryScreen_metadata %>%
  dplyr::select(geneName, DPtier, VemRtier)
primaryScreen_metadata$geneName = gsub("BRD2-BD1", "BRD2", primaryScreen_metadata$geneName)
primaryScreen_metadata$geneName = gsub("BRD2-BD2", "BRD2", primaryScreen_metadata$geneName)
primaryScreen_metadata = primaryScreen_metadata %>%
  distinct()

#Add data to master table 
primaryScreenResults = left_join(primaryScreenResults_DP, primaryScreenResults_VemR, by = c("gRNA", "library"))
masterTable = left_join(masterTable, primaryScreenResults, by = c("gRNA", "library"))
masterTable = left_join(masterTable, primaryScreen_metadata, by = c("gene" = "geneName"))


#Load follow-up screen data
secondaryScreen_targets = read.table("02_metadata/TargetedScreen_targets.txt", header = T)
secondaryScreen_targets$IncludedInSecondaryScreen = "yes"
secondaryScreen_targets$target = NULL

SecondaryScreen_results_WM989_state = read.table("03_extractedData/190122_targetedCrisprScreen_hitSelection/DPanalysis_WM989.bulk_allInfo.txt", header = T)
SecondaryScreen_results_451_state = read.table("03_extractedData/190122_targetedCrisprScreen_hitSelection/DPanalysis_X451Lu_allInfo.txt", header = T)
SecondaryScreen_results_WM989_fate = read.table("03_extractedData/190122_targetedCrisprScreen_hitSelection/VemRanalysis_WM989.bulk_allInfo.txt", header = T)
SecondaryScreen_results_451_fate = read.table("03_extractedData/190122_targetedCrisprScreen_hitSelection/VemRanalysis_X451Lu_allInfo.txt", header = T)

SecondaryScreen_results_WM989_state = SecondaryScreen_results_WM989_state %>%
  dplyr::select(gRNA, lFC, medianlFC, sdlFC)
colnames(SecondaryScreen_results_WM989_state) = c("gRNA", "secondaryScreen_state_WM989-A6-G3-Cas9_lFC", "secondaryScreen_state_WM989-A6-G3-Cas9_medianlFC", "secondaryScreen_state_WM989-A6-G3-Cas9_sdlFC")

SecondaryScreen_results_451_state = SecondaryScreen_results_451_state %>%
  dplyr::select(gRNA, lFC, medianlFC, sdlFC)
colnames(SecondaryScreen_results_451_state) = c("gRNA", "secondaryScreen_state_451Lu-Cas9_lFC", "secondaryScreen_state_451Lu-Cas9_medianlFC", "secondaryScreen_state_451Lu-Cas9_sdlFC")

SecondaryScreen_results_WM989_fate = SecondaryScreen_results_WM989_fate %>%
  dplyr::select(gRNA, lFC, medianlFC, sdlFC)
colnames(SecondaryScreen_results_WM989_fate) = c("gRNA", "secondaryScreen_fate_WM989-A6-G3-Cas9_lFC", "secondaryScreen_fate_WM989-A6-G3-Cas9_medianlFC", "secondaryScreen_fate_WM989-A6-G3-Cas9_sdlFC")

SecondaryScreen_results_451_fate = SecondaryScreen_results_451_fate %>%
  dplyr::select(gRNA, lFC, medianlFC, sdlFC)
colnames(SecondaryScreen_results_451_fate) = c("gRNA", "secondaryScreen_fate_451Lu-Cas9_lFC", "secondaryScreen_fate_451Lu-Cas9_medianlFC", "secondaryScreen_fate_451Lu-Cas9_sdlFC")

secondaryScreen_results = left_join(secondaryScreen_targets, SecondaryScreen_results_WM989_state, by = "gRNA") %>%
  left_join(., SecondaryScreen_results_WM989_fate, by = "gRNA") %>%
  left_join(., SecondaryScreen_results_451_state, by = "gRNA") %>%
  left_join(., SecondaryScreen_results_451_fate, by = "gRNA")

#Add data to master table
masterTable = left_join(masterTable, secondaryScreen_results, by = "gRNA")

#Load IF data
IF_allResults = read.table("03_extractedData/180925_NGFR-IF/allResults.txt", header = T)
IF_allResults = IF_allResults %>%
  dplyr::select(gRNA, numCells, log2_NGFRhi_FC, meanlFC, sdlFC)
colnames(IF_allResults) = c("gRNA", "IF_WM989-A6-G3-Cas9-5a3_numberOfCells", "IF_WM989-A6-G3-Cas9-5a3_FreqOfNGFRhiCells_lFC","IF_WM989-A6-G3-Cas9-5a3_FreqOfNGFRhiCells_meanlFC", "IF_WM989-A6-G3-Cas9-5a3_FreqOfNGFRhiCells_sdlFC")

IF_allResults = IF_allResults %>%
  mutate(IncludedInIF = "yes")

#Add data to master Table
masterTable = left_join(masterTable, IF_allResults, by = "gRNA")

#Load sgRNAs used
colonyGrowth_guides = read.table("02_metadata/181101_ColonyGrowth36_validationOfKOs_reanalyzed/sgRNAsUsed.txt", header = T)
colonyGrowth_results = read.table("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/colonyGrowthResults_allhits.txt", header = T)

colonyGrowth_results = colonyGrowth_results %>%
  dplyr::select(target, Rcolonies_lFC)

colonyGrowth_results = left_join(colonyGrowth_guides, colonyGrowth_results, by = c("gene" = "target"))

colonyGrowth_results = colonyGrowth_results %>%
  dplyr::select(-gene)
colnames(colonyGrowth_results) = c("gRNA", "colonyGrowth_WM989-A6-G3-Cas9-5a3_lFC")

colonyGrowth_results = colonyGrowth_results %>%
  mutate(IncludedInColonyGrowth = "yes")

#Add data to master table
masterTable = left_join(masterTable, colonyGrowth_results, by = "gRNA")

#Load RNA-seq samples
RNAseqSamples = read.table("02_metadata/TargetedScreen_targets.txt", header = T)
RNAseqSamples = RNAseqSamples %>%
  dplyr::select(-target) %>%
  mutate(IncludedInRNAseq = "yes")

#Add data to master table
masterTable = left_join(masterTable, RNAseqSamples, by = "gRNA")

#Erase unnecessary fields
masterTable = masterTable %>%
  dplyr::select(-target)

# save table
write.table(masterTable, "03_extractedData/190430_makeMasterTable/190430_masterTable.txt", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE)

