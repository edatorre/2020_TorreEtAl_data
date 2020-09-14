# #make validation tables based on IF and colony growth
# rm(list = ls())

#Load name conversion file
nameKey = read.table("02_metadata/190402_screenTargets_geneDomainTarget_nameConversion.txt", header = T)
nameKey = nameKey %>%
  dplyr::select(-gRNA) %>%
  distinct()

#Define the positive and negative controls
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
negativeControls = c("WT", "neg")

#Load primary screen metadata
metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)
metadata = metadata %>%
  dplyr::select(geneName, EffectOnScreen_DP, EffectOnScreen_VemR) %>%
  distinct()
metadata = left_join(metadata, nameKey, by = c("geneName" = "target"))

#Load IF summarized results
IF_summarizedResults = read.table("03_extractedData/180925_NGFR-IF/sumarizedResults.txt", header = T)

#Remove controls from summary table
IF_summarizedResults = IF_summarizedResults %>%
  filter(!target %in% positiveControls) %>%
  filter(!target %in% negativeControls) %>%
  distinct()

#Load colony growth data
CG_summarizedResults = read.table("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/colonyGrowthResults_allhits.txt", header = T)

#Select relevant fields
CG_summarizedResults = CG_summarizedResults %>%
  dplyr::select(target, Rcolonies_lFC) %>%
  distinct()

#Combine results from 
IF_results = left_join(metadata, IF_summarizedResults, by = c("geneName" = "target"))
# metadata$geneName = gsub("BRD2-BD1", "BRD2", metadata$geneName)
CGresults = left_join(metadata, CG_summarizedResults, by = c("gene" = "target"))

####Make list of hits that validate
validationTable_DP = IF_results %>%
  mutate(IF_validationStatus = ifelse((EffectOnScreen_DP == "down" & meanlFC <= -0.5) | (EffectOnScreen_DP == "up" & meanlFC >= 0.5), "tested_validated", "tested_NotValidated"))
validationTable_DP$IF_validationStatus[is.na(validationTable_DP$IF_validationStatus)] = "notTested"

validationTable_VemR = CGresults %>%
  mutate(CG_validationStatus = ifelse((EffectOnScreen_VemR == "down" & Rcolonies_lFC <= -0.5) | (EffectOnScreen_VemR == "up" & Rcolonies_lFC >= 0.5), "tested_validated", "tested_NotValidated"))
validationTable_VemR$CG_validationStatus[is.na(validationTable_VemR$CG_validationStatus)] = "notTested"

#save table
write.table(validationTable_DP, file = "03_extractedData/190406_validationTables/validationTable_DP.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(validationTable_VemR, file = "03_extractedData/190406_validationTables/validationTable_VemR.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

