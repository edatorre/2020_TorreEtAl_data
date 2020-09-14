# #The goal of this script is to compare frequency of NGFR high cells to Vemurafenib resistance given a KO.

#Select controls to eliminate from analysis
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
negativeControls = c("neg", "WT", "Neg01", "Neg02")

#Load metadata
metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)
metadata = metadata %>%
  filter(!geneName %in% positiveControls) %>%
  filter(!geneName %in% negativeControls)
metadata$geneName = gsub("BRD2-BD2", "BRD2", metadata$geneName)

#Keep only metadata needed
metadata = metadata %>%
  filter(!(geneName == "EP300" & DPtier == "tier066")) %>%
  filter(!(geneName == "PBRM1" & DPtier == "tierNotHit")) %>%
  filter(!(geneName == "PBRM1" & library == "epigenetic"))

#Load cluster information
clusterGroupings = read.table("03_extractedData/190321_KO-RNAseq/GSEA/clusterContents/c5.bp.v6.2_clusterContents_KOs.txt", header = T)
tmp = clusterGroupings %>%
  dplyr::select(row_dividedTree) %>%
  distinct()

tmp = tmp %>%
  mutate(TtoB = 1:nrow(tmp)) %>%
  mutate(cluster = paste("cluster_", TtoB, sep = ""))

clusterGroupings = left_join(clusterGroupings, tmp, by = "row_dividedTree")
clusterGroupings$row_dividedTree = NULL
clusterGroupings$TtoB = NULL
clusterGroupings$KO = gsub("BRD2.BD2", "BRD2", clusterGroupings$KO)
  
#Load NGFR IF data.
NGFRif = read.table(file = "03_extractedData/180925_NGFR-IF/sumarizedResults.txt", header = T)
colnames(NGFRif)[5] = "meanlFC_IF"
colnames(NGFRif)[7] = "selFC_IF"

#Eliminate one of the BRD2 samples from the NGFR-IF dataser (I kept the one with the bromodomain I targeted in the colony groth assay. It has the weakest effect of the two.)
NGFRif$target = gsub("BRD2-BD2", "BRD2", NGFRif$target)

#Load Resistance data
ResistanceData = read.table( file = "03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/colonyGrowthResults_allhits.txt", header = T)
ResistanceData = ResistanceData %>%
  dplyr::select(target, Rcolonies_lFC)

#merge IF and colony growth data
mergedData = left_join(ResistanceData, NGFRif, by = "target")

#add metadata
mergedData = left_join(mergedData, metadata, by = c("target" = "geneName"))

#Add cluster groupings
mergedData = left_join(mergedData, clusterGroupings, by = c("target" = "KO"))

#Add screen of origin
#Highlight which screen the target came from
  #Convert tiers into numbers to use later
  mergedData$DPtier = gsub("tier075", 1, mergedData$DPtier)
  mergedData$DPtier = gsub("tier066", 2, mergedData$DPtier)
  mergedData$DPtier = gsub("tier050", 3, mergedData$DPtier)
  mergedData$DPtier = gsub("tierNotHit", 4, mergedData$DPtier)
  
  mergedData$VemRtier = gsub("tier075", 1, mergedData$VemRtier)
  mergedData$VemRtier = gsub("tier066", 2, mergedData$VemRtier)
  mergedData$VemRtier = gsub("tier050", 3, mergedData$VemRtier)
  mergedData$VemRtier = gsub("tierNotHit", 4, mergedData$VemRtier)
  
  #Detemine color of label
  mergedData$labelColor = "none"
  
  #Determine colors for labels (I select the color based on the the screen where the target is either tier 1 or tier 2. If that target is tier 1 or 2 in both screens then I call it a hit in both screens)
  mergedData$labelColor = ifelse((mergedData$DPtier == 1 | mergedData$DPtier == 2), "stateHit", mergedData$labelColor)
  mergedData$labelColor = ifelse((mergedData$VemRtier == 1 | mergedData$VemRtier == 2), "fateHit", mergedData$labelColor)
  mergedData$labelColor = ifelse((mergedData$DPtier == 1 | mergedData$DPtier == 2) & (mergedData$VemRtier == 1 | mergedData$VemRtier == 2), "hitInBothScreen", mergedData$labelColor)
  mergedData$labelColor = ifelse((mergedData$DPtier == 3 | mergedData$DPtier == 4) & (mergedData$VemRtier == 3 | mergedData$VemRtier == 4), "tier3And4Targets", mergedData$labelColor)

#Remove controls
mergedData = mergedData %>%
  filter(!target %in% negativeControls) %>%
  filter(!target %in% positiveControls)

#make tmp table with swithout the individual sgRNA lFC, otherwise the KO names appear multiple times
tmp = mergedData %>%
  dplyr::select(-gRNA, -log2_NGFRhi_FC) %>%
  distinct()

#Make plots
plot_IFoverColonies = ggplot(data = tmp, aes(x = meanlFC_IF, y = Rcolonies_lFC, color = labelColor, size = 2)) +
  geom_point() +
  #geom_text(aes(label = target), nudge_x = 0.25, nudge_y = 0.25) +
  geom_text_repel(aes(label = target), nudge_x = 0.25, nudge_y = 0.25) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  #geom_smooth(method = "lm") +
  theme_classic()
plot_IFoverColonies
ggsave("04_plots/190218_NGFR-IFvsResistance/IFoverColonies_byScreenOfOrigin.PDF")

plot_IFoverColonies = ggplot(data = tmp, aes(x = meanlFC_IF, y = Rcolonies_lFC, color = cluster, size = 2)) +
  geom_point() +
  geom_text_repel(aes(label = target), nudge_x = 0.25, nudge_y = 0.25) +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  #geom_smooth(method = "lm") +
  theme_classic()
plot_IFoverColonies
ggsave("04_plots/190218_NGFR-IFvsResistance/IFoverColonies_byCluster.PDF")
