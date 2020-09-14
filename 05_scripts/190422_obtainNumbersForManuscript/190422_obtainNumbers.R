#This script aims to obtain all the numbers used throuhgout the manuscript. 

#Create file to store numbers
file.create("03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt")

#Name controls
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
negativeControls = c("neg", "WT")

#Load name conversion file
nameKey = read.table("02_metadata/screenTargets_geneDomainTarget_nameConversion.txt", header = T)
nameKey = nameKey %>%
  dplyr::select(domain, target) %>%
  distinct()

#### Primary screen targets ####
commentary = "***Primary Screen Targets***"
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#sgRNAs and targets used in the primary screen
primaryScreen_targets = read.table("02_metadata/ScreenTargets_metadata.txt", header = T)
primaryScreen_targets_ControlsRemoved = primaryScreen_targets %>%
  filter(!domain %in% negativeControls) %>%
  filter(!domain %in% positiveControls)

primaryScreen_targets_ControlsRemoved = primaryScreen_targets_ControlsRemoved %>%
  group_by(domain) %>%
  mutate(nGuidesPerDomain = n_distinct(gRNA)) %>%
  ungroup(domain)

primaryScreen_targets_ControlsRemoved = primaryScreen_targets_ControlsRemoved %>%
  filter(nGuidesPerDomain > 3)

primaryScreen_targets_domains = primaryScreen_targets_ControlsRemoved %>%
  dplyr::select(gene, domain, library) %>%
  distinct()

primaryScreen_targets_genes = primaryScreen_targets_domains %>%
  dplyr::select(gene, library) %>%
  distinct()

primaryScreen_targets_domains_epigenetic = primaryScreen_targets_domains %>%
  filter(library == "epigenetic") %>%
  distinct()
primaryScreen_targets_domains_kinase = primaryScreen_targets_domains %>%
  filter(library == "kinase") %>%
  distinct()
primaryScreen_targets_domains_TF = primaryScreen_targets_domains %>%
  filter(library == "TF") %>%
  distinct()

primaryScreen_targets_genes_epigenetic = primaryScreen_targets_genes %>%
  filter(library == "epigenetic") %>%
  distinct()
primaryScreen_targets_genes_kinase = primaryScreen_targets_genes %>%
  filter(library == "kinase") %>%
  distinct()
primaryScreen_targets_genes_TF = primaryScreen_targets_genes %>%
  filter(library == "TF") %>%
  distinct()

  
commentary = paste("There are", nrow(primaryScreen_targets_ControlsRemoved), "sgRNAs in the libraries (without controls)")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("There are", nrow(primaryScreen_targets_domains), "domain targets in the libraries")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("There are", nrow(primaryScreen_targets_genes), "gene targets in the libraries")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("There are", nrow(primaryScreen_targets_domains_epigenetic), "epigenetic domain targets aimed at", nrow(primaryScreen_targets_genes_epigenetic), "different proteins")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("There are", nrow(primaryScreen_targets_domains_kinase), "kinase domain targets aimed at", nrow(primaryScreen_targets_genes_kinase), "different proteins")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("There are", nrow(primaryScreen_targets_domains_TF), "TF domain targets aimed at", nrow(primaryScreen_targets_genes_TF), "different proteins")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)


#### Primary Screen Results####
commentary = "***Primary Screen Results***"
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#Number of hits in state screen, and directionality of each
DPscreen_epigeneticHits = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_epigenetic_hitList.txt", header = T)
DPscreen_kinaseHits = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_kinase_hitList.txt", header = T)
DPscreen_TFHits = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_TF_hitList.txt", header = T)
DPscreen_hits = bind_rows(DPscreen_epigeneticHits, DPscreen_kinaseHits, DPscreen_TFHits)

DPscreen_hits = DPscreen_hits %>%
  filter(!domain %in% negativeControls) %>%
  filter(!domain %in% positiveControls)

DPscreen_hits_up = DPscreen_hits %>%
  filter(medianlFC > 0)
DPscreen_hits_down = DPscreen_hits %>%
  filter(medianlFC < 0)

commentary = paste("There are", nrow(DPscreen_hits), "hits in the primary state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste(nrow(DPscreen_hits_up), "of the hits increase the frequency of NGFRhi EGFRhi cells")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste(nrow(DPscreen_hits_down), "of the hits decrease the frequency of NGFRhi EGFRhi cells")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#Number of hits in fate screen, and directionality of each
VemRscreen_epigeneticHits = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_epigenetic_hitList.txt", header = T)
VemRscreen_kinaseHits = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_kinase_hitList.txt", header = T)
VemRscreen_TFHits = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_TF_hitList.txt", header = T)
VemRscreen_hits = bind_rows(VemRscreen_epigeneticHits, VemRscreen_kinaseHits, VemRscreen_TFHits)

VemRscreen_hits_up = VemRscreen_hits %>%
  filter(medianlFC > 0)
VemRscreen_hits_down = VemRscreen_hits %>%
  filter(medianlFC < 0)

commentary = paste("There are", nrow(VemRscreen_hits), "hits in the primary fate screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste(nrow(VemRscreen_hits_up), "of the hits increase the frequency of cells resistant to vemurafenib")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste(nrow(VemRscreen_hits_down), "of the hits decrease the frequency of cells resistant to vemurafenib")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#Overlap between screens
overlapBetweenPrimaryScreens = inner_join(DPscreen_hits, VemRscreen_hits, by = "domain")
  
commentary = paste(nrow(overlapBetweenPrimaryScreens), "targets are hits in both state and fate screens")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)


#### Targeted Screen Composition ####
commentary = "***Secondary Screen Targets***"
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#Number of targets on targted screen
targetedScreen_sgRNAs = read.table("02_metadata/TargetedScreen_targets.txt", header = T)
targetedScreen_sgRNAs = left_join(targetedScreen_sgRNAs, nameKey, by = "target")
targetedScreen_targets = targetedScreen_sgRNAs %>%
  dplyr::select(target, domain) %>%
  filter(!target %in% positiveControls | !target %in% negativeControls) %>%
  distinct()

#Number of controls included
targetedScreen_controls = targetedScreen_sgRNAs %>%
  filter(target %in% positiveControls | target %in% negativeControls)

#Number of hits from each screen used in targeted screen
targetedScreen_fromDPhits = DPscreen_hits %>%
  filter(domain %in% targetedScreen_targets$domain)
targetedScreen_fromVemRhits = VemRscreen_hits %>%
  filter(domain %in% targetedScreen_targets$domain)

targetedScreen_inlcudedHitInBothScreens = inner_join(targetedScreen_fromDPhits, targetedScreen_fromVemRhits, by = "domain")

#total targets included because they are tier 1 or tier 2
targetedScreen_hitsFromPrimaryScreens = bind_rows(targetedScreen_fromDPhits, targetedScreen_fromVemRhits) 
targetedScreen_hitsFromPrimaryScreens = targetedScreen_hitsFromPrimaryScreens %>%
  dplyr::select(domain) %>%
  distinct()

#targets in tiers 3 and tiers 4
targetedScreen_notHits = targetedScreen_targets %>%
  filter(!domain %in% targetedScreen_hitsFromPrimaryScreens$domain) %>%
  distinct()

tmp = inner_join(targetedScreen_hitsFromPrimaryScreens, targetedScreen_notHits, by = "domain")

commentary = paste("There are",paste(nrow(targetedScreen_sgRNAs), "sgRNAs in the targeted screen aimed at", nrow(targetedScreen_targets),"different targets"))
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("The targeted screen includes", nrow(targetedScreen_controls), "control sgRNAs")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste(nrow(targetedScreen_hitsFromPrimaryScreens), "of the targets are tier 1 or tier 2 targets in at least one of the primary screens")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("The other",nrow(targetedScreen_notHits), "targets are tier 3 or tier 4 targets in both screens")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste(paste(nrow(targetedScreen_fromDPhits), "of the targets are hits in the state screen"))
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste(paste(nrow(targetedScreen_fromVemRhits), "of the targets are hits in the fate screen"))
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste(paste(nrow(targetedScreen_inlcudedHitInBothScreens), "of the targets above are hits in both screens"))
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#### Targeted Screen Results ####
commentary = "***Secondary Screen Results***"
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)
  
#Number of hits validated in the targeted screen: DP on WM989
targetedScreen_results_DP_WM989 = read.table("03_extractedData/190122_targetedCrisprScreen_hitSelection/DPanalysis_WM989.bulk_hitList.txt", header = T)

tmp = DPscreen_hits %>%
  dplyr::select(domain, medianlFC) %>%
  distinct()
colnames(tmp) = c("domain", "primaryScreen_medianlFC")

targetedScreen_results_DP_WM989 = left_join(targetedScreen_results_DP_WM989, nameKey, by = "target")
targetedScreen_results_DP_validatedFromPrimary_WM989 = inner_join(targetedScreen_results_DP_WM989, tmp, by = "domain")
targetedScreen_results_DP_validatedFromPrimary_WM989 = targetedScreen_results_DP_validatedFromPrimary_WM989 %>%
  filter((primaryScreen_medianlFC > 0 & medianlFC > 0) | (primaryScreen_medianlFC < 0 & medianlFC < 0)) %>%
  filter(!domain %in% negativeControls | !domain %in% positiveControls)

commentary = paste(nrow(targetedScreen_results_DP_validatedFromPrimary_WM989), "of the", nrow(targetedScreen_fromDPhits) ," state Tier 1 and Tier 2 targets from the primary state screen also show a two-fold change in the WM989 cells in the targeted screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#Number of hits validated in the targeted screen: DP on 451Lu
targetedScreen_results_DP_451lu = read.table("03_extractedData/190122_targetedCrisprScreen_hitSelection/DPanalysis_X451Lu_hitList.txt", header = T)

tmp = DPscreen_hits %>%
  dplyr::select(domain, medianlFC) %>%
  distinct()
colnames(tmp) = c("domain", "primaryScreen_medianlFC")

targetedScreen_results_DP_451Lu = left_join(targetedScreen_results_DP_451lu, nameKey, by = "target")
targetedScreen_results_DP_validatedFromPrimary_451Lu = inner_join(targetedScreen_results_DP_451Lu, tmp, by = "domain")
targetedScreen_results_DP_validatedFromPrimary_451Lu = targetedScreen_results_DP_validatedFromPrimary_451Lu %>%
  filter((primaryScreen_medianlFC > 0 & medianlFC > 0) | (primaryScreen_medianlFC < 0 & medianlFC < 0)) %>%
  filter(!domain %in% negativeControls | !domain %in% positiveControls)

commentary = paste(nrow(targetedScreen_results_DP_validatedFromPrimary_451Lu), "of the", nrow(targetedScreen_fromDPhits) ," state Tier 1 and Tier 2 targets from the primary state screen also show a two-fold change in the 451Lu cells in the targeted screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#Number of tagret in the targeted screen that replicated in both cell lines
targetedScreen_DP_doubleRepilation = targetedScreen_results_DP_validatedFromPrimary_451Lu %>%
  filter(domain %in% targetedScreen_results_DP_validatedFromPrimary_WM989$domain)

#Number of hits validated in the targeted screen: VemR on WM989
targetedScreen_results_VemR_WM989 = read.table("03_extractedData/190122_targetedCrisprScreen_hitSelection/VemRanalysis_WM989.bulk_allInfo.txt", header = T)
targetedScreen_results_VemR_WM989 = targetedScreen_results_VemR_WM989 %>%
  dplyr::select(target, medianlFC, nGuidesOverFC) %>%
  distinct()

tmp = VemRscreen_hits %>%
  dplyr::select(domain, medianlFC) %>%
  distinct()
colnames(tmp) = c("domain", "primaryScreen_medianlFC")

targetedScreen_results_VemR_WM989 = left_join(targetedScreen_results_VemR_WM989, nameKey, by = "target")
targetedScreen_results_VemR_validatedFromPrimary_WM989 = inner_join(targetedScreen_results_VemR_WM989, tmp, by = "domain")
targetedScreen_results_VemR_validatedFromPrimary_WM989 = targetedScreen_results_VemR_validatedFromPrimary_WM989 %>%
  filter((primaryScreen_medianlFC > 0 & medianlFC > 0) | (primaryScreen_medianlFC < 0 & medianlFC < 0)) %>%
  filter(!domain %in% negativeControls | !domain %in% positiveControls)

commentary = paste(nrow(targetedScreen_results_VemR_validatedFromPrimary_WM989), "of the", nrow(targetedScreen_fromVemRhits) ," fate Tier 1 and Tier 2 targets from the primary state screen also show a two-fold change in the WM989 cells in the targeted screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#Number of hits validated in the targeted screen: VemR on 451Lu
targetedScreen_results_VemR_451Lu = read.table("03_extractedData/190122_targetedCrisprScreen_hitSelection/VemRanalysis_X451Lu_allInfo.txt", header = T)
targetedScreen_results_VemR_451Lu = targetedScreen_results_VemR_451Lu %>%
  dplyr::select(target, medianlFC, nGuidesOverFC) %>%
  distinct()

tmp = VemRscreen_hits %>%
  dplyr::select(domain, medianlFC) %>%
  distinct()
colnames(tmp) = c("domain", "primaryScreen_medianlFC")

targetedScreen_results_VemR_451Lu = left_join(targetedScreen_results_VemR_451Lu, nameKey, by = "target")
targetedScreen_results_VemR_validatedFromPrimary_451Lu = inner_join(targetedScreen_results_VemR_451Lu, tmp, by = "domain")
targetedScreen_results_VemR_validatedFromPrimary_451Lu = targetedScreen_results_VemR_validatedFromPrimary_451Lu %>%
  filter((primaryScreen_medianlFC > 0 & medianlFC > 0) | (primaryScreen_medianlFC < 0 & medianlFC < 0)) %>%
  filter(!domain %in% negativeControls | !domain %in% positiveControls)

commentary = paste(nrow(targetedScreen_results_VemR_validatedFromPrimary_451Lu), "of the", nrow(targetedScreen_fromVemRhits) ," fate Tier 1 and Tier 2 targets from the primary state screen also show a two-fold change in the 451Lu cells in the targeted screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#Number of tagret in the targeted screen that replicated in both cell lines
targetedScreen_VemR_doubleRepilation = targetedScreen_results_VemR_validatedFromPrimary_451Lu %>%
  filter(domain %in% targetedScreen_results_VemR_validatedFromPrimary_WM989$domain)

#### IF validation ####
commentary = "***IF validation***"
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#Number of targtes from each screen validated by IF

IF_metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)
IF_metadata = IF_metadata %>%
  dplyr::select(geneName, medianlFC_DP, medianlFC_VemR, DPtier, VemRtier)
colnames(IF_metadata) = c("geneName", "primaryScreen_medianlFC_DP", "primaryScreen_medianlFC_VemR", "DPtier", "VemRtier")

IF_metadata = left_join(nameKey, IF_metadata, by = c("target" = "geneName"))

# IF_plateMap = read.table("02_metadata/180925_NGFR-IF/NGFR-IF_plateMaps.txt", header = T)
IF_plateMap = read.table("02_metadata/180925_NGFR-IF/NGFR-IF_sampleKey.txt", header = T)
IF_targets = IF_plateMap %>%
  dplyr::select(target) %>%
  distinct()
IF_results = read.table("03_extractedData/180925_NGFR-IF/sumarizedResults_targetsOnly.txt", header = T)
IF_controls = IF_plateMap %>%
  filter(target %in% negativeControls | target %in% positiveControls)

IF_results = left_join(IF_metadata, IF_results, by = "target")

#****** I think the problem is that the IF results have too many rows coming in...****

#Rmeove duplicated data
IF_results = IF_results %>%
  filter(!(target == "PBRM1" & VemRtier == "tierNotHit")) %>%
  filter(!(target == "EP300" & DPtier == "tier066"))

commentary = paste("We assayed the frequency of NGFR-high cells in", nrow(IF_plateMap), "samples (including controls), targeting", nrow(IF_targets) ,"different proteins")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("The dataset included", nrow(IF_controls), "controls")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

# Number of targets from tiers one and two of state screen images
IF_results_tier1and2 = IF_results %>%
  filter(DPtier == "tier066" | DPtier == "tier075")
IF_results_tier3and4 = IF_results %>%
  filter(DPtier == "tier050" | DPtier == "tierNotHit")

commentary = paste("Of those",nrow(IF_targets), "targets,", nrow(IF_results_tier1and2), "were tier one or tier 2 targets in the state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#Number of targets (by tier) showing a 50% change in frequency of NGFRhi cells concordant with the primary screen
IF_results_tier1and2_validated_50percentChange = IF_results_tier1and2 %>%
  filter((primaryScreen_medianlFC_DP > 0 & meanlFC >= 0.5) | (primaryScreen_medianlFC_DP < 0 & meanlFC <= -0.5))
IF_results_tier3and4_validated_50percentChange = IF_results_tier3and4 %>%
  filter((primaryScreen_medianlFC_DP > 0 & meanlFC >= 0.5) | (primaryScreen_medianlFC_DP < 0 & meanlFC <= -0.5))

IF_results_tier1and2_validated_twoFoldChange = IF_results_tier1and2 %>%
  filter((primaryScreen_medianlFC_DP > 0 & meanlFC >= 1) | (primaryScreen_medianlFC_DP < 0 & meanlFC <= -1))
IF_results_tier3and4_validated_twoFoldChange = IF_results_tier3and4 %>%
  filter((primaryScreen_medianlFC_DP > 0 & meanlFC >= 1) | (primaryScreen_medianlFC_DP < 0 & meanlFC <= -1))

commentary = paste("Of the",nrow(IF_results_tier1and2), "tier1 and tier 2 targets,", nrow(IF_results_tier1and2_validated_50percentChange), "showed a 50% change in frequency of NGFRhi cells conconrdant with the primary state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("Of the",nrow(IF_results_tier1and2), "tier1 and tier 2 targets,", nrow(IF_results_tier1and2_validated_twoFoldChange), "showed a two-fold change in frequency of NGFRhi cells conconrdant with the primary state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("Of the",nrow(IF_results_tier3and4), "tier 3 and tier 4 targets,", nrow(IF_results_tier3and4_validated_50percentChange), "showed a 50% change in frequency of NGFRhi cells conconrdant with the primary state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("Of the",nrow(IF_results_tier3and4), "tier3 and tier 4 targets,", nrow(IF_results_tier3and4_validated_twoFoldChange), "showed a 50% change in frequency of NGFRhi cells conconrdant with the primary state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

#### Colony growth validation ####
commentary = "***Colony growth validation***"
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

CG_data = read.table("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/colonyGrowthResults_allhits.txt", header = T)

commentary = paste("We quantified the number of colonies growing under vemurafenib in", nrow(CG_data), "samples (including two non-targeting controls)")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

CG_stateHits = CG_data %>%
  filter(DPtier == "tier1" | DPtier == "tier2")

CG_fateHits = CG_data %>%
  filter(VemRtier == "tier1" | DPtier == "tier2")

CG_stateHits_validated_50percentChange = CG_stateHits %>%
  filter((EffectOnScreen_DP == "up" & Rcolonies_lFC >= 0.5) | (EffectOnScreen_DP == "down" & Rcolonies_lFC <= -0.5))

CG_stateHits_validated_twoFoldChange = CG_stateHits %>%
  filter((EffectOnScreen_DP == "up" & Rcolonies_lFC >= 1) | (EffectOnScreen_DP == "down" & Rcolonies_lFC <= -1))

CG_fateHits_validated_50percentChange = CG_fateHits %>%
  filter((EffectOnScreen_VemR == "up" & Rcolonies_lFC >= 0.5) | (EffectOnScreen_VemR == "down" & Rcolonies_lFC <= -0.5))

CG_fateHits_validated_twoFoldChange = CG_fateHits %>%
  filter((EffectOnScreen_VemR == "up" & Rcolonies_lFC >= 1) | (EffectOnScreen_VemR == "down" & Rcolonies_lFC <= -1))

commentary = paste("Of the", nrow(CG_data), "samples,", nrow(CG_stateHits), "were tier 1 or tier 2 targets from the state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("Of those", nrow(CG_stateHits), "samples,", nrow(CG_stateHits_validated_50percentChange), "showed at least a 50% change in the number of resistant colonies concordant in direction with the primary state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("Of the", nrow(CG_data), "samples,", nrow(CG_fateHits), "were tier 1 or tier 2 targets from the state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("Of those", nrow(CG_fateHits), "samples,", nrow(CG_fateHits_validated_50percentChange), "showed at least a 50% change in the number of resistant colonies concordant in direction with the primary state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)


#### RNA-seq ####
commentary = "***RNA-seq***"
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

RNAseq_metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)

RNAseq_metadata = RNAseq_metadata %>%
  dplyr::select(geneName, medianlFC_DP, medianlFC_VemR, DPtier, VemRtier)
colnames(RNAseq_metadata) = c("geneName", "primaryScreen_medianlFC_DP", "primaryScreen_medianlFC_VemR", "DPtier", "VemRtier")

RNAseqGuides = read.table("02_metadata/TargetedScreen_targets.txt", header = T)
RNAseqControls = RNAseqGuides %>%
  filter(target == "neg")
RNAseqTargets = RNAseqGuides %>%
  dplyr::select(target) %>%
  filter(!target %in% negativeControls) %>%
  #filter(!target %in% positiveControls) %>%
  distinct()

RNAseqTargets = left_join(RNAseqTargets, RNAseq_metadata, by = c("target" = "geneName"))

#Rmeove duplicated data
RNAseqTargets = RNAseqTargets %>%
  filter(!(target == "PBRM1" & VemRtier == "tierNotHit")) %>%
  filter(!(target == "EP300" & DPtier == "tier066"))

# Targets in tier 1 and tier 2 of state screen
RNAseqTargets_stateHits = RNAseqTargets %>%
  filter(DPtier == "tier075" | DPtier == "tier066")

# Targets in tier 1 and tier 2 of fate screen
RNAseqTargets_fateHits = RNAseqTargets %>%
  filter(VemRtier == "tier075" | VemRtier == "tier066")

#Hits in both screens
RNAseqTargets_fromBothScreens = RNAseqTargets_stateHits %>%
  filter(target %in% RNAseqTargets_fateHits$target)

# Targets not in tier 1 or tier 2 in any of the screens
RNAseqTargets_notHits = RNAseqTargets %>%
  filter(!target %in% RNAseqTargets_stateHits$target) %>%
  filter(!target %in% RNAseqTargets_fateHits$target)


commentary = paste("We assayed the trasncriptome of", nrow(RNAseqGuides), "samples, including", nrow(RNAseqControls), "non-targeting controls")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("In total, we targeted", nrow(RNAseqTargets), "different proteins")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("Of the", nrow(RNAseqTargets), "targets,", nrow(RNAseqTargets_stateHits) , "are tier 1 or tier 2 targets in the state screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("Of the", nrow(RNAseqTargets), "targets,", nrow(RNAseqTargets_fateHits) , "are tier 1 or tier 2 targets in the fate screen")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("Of the", nrow(RNAseqTargets), "targets,", nrow(RNAseqTargets_fromBothScreens) , "are a hit in both screens")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

commentary = paste("Of the", nrow(RNAseqTargets), "targets,", nrow(RNAseqTargets_notHits) , "are in tiers 3 and 4 in both screens")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)

# Number of this that validated in one way or another, regardless of tier. 
targetsValidatedInIF = IF_results %>%
  filter((primaryScreen_medianlFC_DP >= 0 & meanlFC >= 0.5) | (primaryScreen_medianlFC_DP <= 0 & meanlFC <= -0.5)) %>%
  dplyr::select(target) %>%
  distinct()

targetsValidatedInCG = CG_data %>%
  filter((EffectOnScreen_DP == "up" & Rcolonies_lFC >= 0.5) | (EffectOnScreen_DP == "down" & Rcolonies_lFC <= -0.5) |
           (EffectOnScreen_VemR == "up" & Rcolonies_lFC >= 0.5) | (EffectOnScreen_VemR == "down" & Rcolonies_lFC <= -0.5)) %>%
  dplyr::select(target) %>%
  distinct()
targetsValidatedInCG$target = gsub("BRD2", "BRD2-BD1", targetsValidatedInCG$target)

targetsValidated = bind_rows(targetsValidatedInIF, targetsValidatedInCG)
targetsValidated = targetsValidated %>%
  distinct()

commentary = paste("Of the", nrow(RNAseqTargets), "targets,", nrow(targetsValidated) , "validated via IF or via colony formation")
print(commentary)
write(commentary, file="03_extractedData/190422_obtainNumbersForManuscript/manuscriptNumbers.txt", append=TRUE)











