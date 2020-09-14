#This script aims to plot the enrichment score of a given KO for a specific gene set vs the number of colonies growth in vemurafenib.

#Load name conversion file
nameConversionFile = read.table("02_metadata/190402_screenTargets_geneDomainTarget_nameConversion.txt", header = T)

#Load the enrichment scores of KOs
enrichmentScores = read.table("03_extractedData/190321_KO-RNAseq/GSEA/ESfiles/c5.bp.v6.2_combinedESfile.txt", header = T)

#Load the colony formation numbers (and fix the BRD2 name to match the sgRNA used on colony formation assay)
resistantColonies = read.table("03_extractedData/181101_ColonyGrowth36_validationOfKOs_reanalyzed/colonyGrowthResults_allhits.txt", header = T)
resistantColonies$target = gsub("BRD2", "BRD2.BD1", resistantColonies$target)

#Load cluster groupings
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

#Keep only KOs that have at least a two-fold change in the number of resistant colonies
resistantColonies = resistantColonies %>%
  filter(Rcolonies_lFC >= 1 | Rcolonies_lFC <= -1)

#Restrict table to relevant fields
resistantColonies = resistantColonies %>%
  dplyr::select(target, Rcolonies_lFC)

#Select KOs
listOfKO = c("SOX10", "MITF", "BRD2.BD1", "DOT1L", "LATS2", "CSK")

#Select gene set 
geneSetList = c("GO_REGULATION_OF_ENDOTHELIAL_CELL_DIFFERENTIATION", "GO_NEGATIVE_REGULATION_OF_CELL_DIFFERENTIATION", "GO_NEURAL_CREST_CELL_DIFFERENTIATION", 
                "GO_STEM_CELL_DIFFERENTIATION", "GO_MESENCHYMAL_CELL_DIFFERENTIATION", "GO_MELANOCYTE_DIFFERENTIATION", "GO_EPITHELIAL_CELL_DIFFERENTIATION", 
                "GO_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION", "GO_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION")
geneSet = "GO_MELANOCYTE_DIFFERENTIATION"

#Limit EnrichmentScore table to KO and gene set of interest
enrichmentScores = enrichmentScores %>%
  filter(GS %in% geneSetList)
enrichmentScores_filtered = enrichmentScores %>%
  #dplyr::select(one_of(listOfKO)) %>%
  filter(GS == "GO_NEURAL_CREST_CELL_DIFFERENTIATION")

#Make plotting table
plottingTable = data.frame(t(enrichmentScores_filtered))
plottingTable$target = rownames(plottingTable)
plottingTable = inner_join(plottingTable, resistantColonies, by = "target")
plottingTable = inner_join(plottingTable, clusterGroupings, by = c("target" = "KO"))
colnames(plottingTable) = c("EnrichmentScoreOnGeneSet", "target", "NumberOfResistantColonies", "cluster")
plottingTable$EnrichmentScoreOnGeneSet = as.numeric(as.character(plottingTable$EnrichmentScoreOnGeneSet))

#Make the plot
resistanceVsDifferentiation = ggplot(data = plottingTable, aes(x = EnrichmentScoreOnGeneSet, y = NumberOfResistantColonies, color = cluster)) +
  geom_point() +
  geom_text_repel(aes(label = target)) +
  ylim(-5, 10) +
  xlim(-2.5, 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Number of colonies resistant to vemurafenib") +
  xlab("Enrichment Score in neural crest differentiation gene set") +
  theme_classic()
resistanceVsDifferentiation
ggsave("04_plots/190321_KO-RNAseq/RessitanceVsDifferentiationStatus.pdf")
