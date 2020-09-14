# #The goal of this script is to do a PCA anlysis on a select number of genes

#Load name conversion key
nameConversionKey = read.table("02_metadata/190321_KO-RNAseq/hg19_nameConvertionFile.txt", header = T)

#Select controls to eliminate from analysis
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
negativeControls = c("neg", "WT")

#Load metadata
metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)
metadata = metadata %>%
  filter(!geneName %in% positiveControls) %>%
  filter(!geneName %in% negativeControls)
metadata$geneName = gsub("-", "\\.", metadata$geneName)

#Keep only metadata needed
metadata = metadata %>%
  filter(!(geneName == "EP300" & DPtier == "tier066")) %>%
  filter(!(geneName == "PBRM1" & DPtier == "tierNotHit")) %>%
  filter(!(geneName == "PBRM1" & library == "epigenetic"))


#Convert tiers into numbers to use later
metadata$DPtier = gsub("tier075", 1, metadata$DPtier)
metadata$DPtier = gsub("tier066", 2, metadata$DPtier)
metadata$DPtier = gsub("tier050", 3, metadata$DPtier)
metadata$DPtier = gsub("tierNotHit", 4, metadata$DPtier)

metadata$VemRtier = gsub("tier075", 1, metadata$VemRtier)
metadata$VemRtier = gsub("tier066", 2, metadata$VemRtier)
metadata$VemRtier = gsub("tier050", 3, metadata$VemRtier)
metadata$VemRtier = gsub("tierNotHit", 4, metadata$VemRtier)

#Detemine color of label
metadata$labelColor = "none"

#Determine colors for labels (I select the color based on the effect on the screen where the targer was the highest tier. If there is a tie, I always select attribute the targets to the screen where it showed the largest effect)
#Select the effect of the higest tier
metadata$labelColor = ifelse(metadata$DPtier <= metadata$VemRtier & metadata$EffectOnScreen_DP == "up", "green_state", metadata$labelColor)
metadata$labelColor = ifelse(metadata$DPtier <= metadata$VemRtier & metadata$EffectOnScreen_DP == "down", "gray_state", metadata$labelColor)
metadata$labelColor = ifelse(metadata$VemRtier < metadata$DPtier & metadata$EffectOnScreen_VemR == "up", "red_fate", metadata$labelColor)
metadata$labelColor = ifelse(metadata$VemRtier < metadata$DPtier & metadata$EffectOnScreen_VemR == "down", "gray_fate", metadata$labelColor)
#When tehre is a tie in tiers, select the strongest effect. 
metadata$labelColor = ifelse((metadata$DPtier == metadata$VemRtier) & (abs(metadata$medianlFC_DP) > abs(metadata$medianlFC_VemR)) &  metadata$medianlFC_DP > 0, "green_state", metadata$labelColor)
metadata$labelColor = ifelse((metadata$DPtier == metadata$VemRtier) & (abs(metadata$medianlFC_DP) > abs(metadata$medianlFC_VemR)) &  metadata$medianlFC_DP < 0, "gray_state", metadata$labelColor)
metadata$labelColor = ifelse((metadata$DPtier == metadata$VemRtier) & (abs(metadata$medianlFC_DP) < abs(metadata$medianlFC_VemR)) &  metadata$medianlFC_VemR > 0, "red_fate", metadata$labelColor)
metadata$labelColor = ifelse((metadata$DPtier == metadata$VemRtier) & (abs(metadata$medianlFC_DP) < abs(metadata$medianlFC_VemR)) &  metadata$medianlFC_VemR < 0, "gray_fate", metadata$labelColor)

#Load enrichment score matrix
ESmatrix = read.table("03_extractedData/190321_KO-RNAseq/GSEA/ESfiles/c5.bp.v6.2_combinedESfile.txt", header = TRUE)

#Eliminate gene symbol column and add rownames
rownames(ESmatrix) = ESmatrix$GS
ESmatrix$GS = NULL

ESmatrix = ESmatrix[rowSums(abs(ESmatrix[1:ncol(ESmatrix)])) > 0 ,]

#Eliminate unwated samples from ESmatrix
ESmatrix = ESmatrix %>%
  dplyr::select(one_of(metadata$geneName))

#Add back the row names
ESmatrix$GS = rownames(ESmatrix)

#Read in enrichment score of new samples: NGFRhi EGFRhi and VemR
enrichmentScoreMatrix_NGFRandEGFRhi = read.table("03_extractedData/190408_RNAseq_otherSamples/NGFRandEGFRhi/GSEA/ESfiles/c5.bp.v6.2_combinedESfile.txt", header = T)
enrichmentScoreMatrix_NGFRandEGFRhi$GS = as.character(enrichmentScoreMatrix_NGFRandEGFRhi$GS)
enrichmentScoreMatrix_otherSamples = read.table("03_extractedData/190408_RNAseq_otherSamples/otherDatatsets/GSEA/ESfiles/c5.bp.v6.2_combinedESfile.txt", header = T)
enrichmentScoreMatrix_otherSamples$GS = as.character(enrichmentScoreMatrix_otherSamples$GS)

#Merge the two new files into a single one
enrichmentScoreMatrix_newSamples = left_join(enrichmentScoreMatrix_NGFRandEGFRhi, enrichmentScoreMatrix_otherSamples, by = "GS")

#Limit the new matrix to only terms in the main matrix
enrichmentScoreMatrix_newSamples = enrichmentScoreMatrix_newSamples %>%
  filter(GS %in% ESmatrix$GS)

#Eliminate gene symbol column and add rownames
rownames(ESmatrix) = ESmatrix$GS
ESmatrix$GS = NULL

rownames(enrichmentScoreMatrix_newSamples) = enrichmentScoreMatrix_newSamples$GS
enrichmentScoreMatrix_newSamples$GS = NULL


#Transpose 
ESmatrix = t(ESmatrix)
enrichmentScoreMatrix_newSamples = t(enrichmentScoreMatrix_newSamples)

#Do PCA
PCAresults = prcomp(ESmatrix, scale. = TRUE)

#Project new samples into PCA
PCAresults_sortedSamples = scale(enrichmentScoreMatrix_newSamples, PCAresults$center, PCAresults$scale) %*% PCAresults$rotation 
PCAresults_sortedSamples = as.data.frame(PCAresults_sortedSamples)

#Extract x
t03_XfromPCA = as.data.frame(PCAresults$x)

#remove row names and add field for the names
t03_XfromPCA$sampleName = rownames(t03_XfromPCA)
rownames(t03_XfromPCA) = NULL
t03_XfromPCA$sampleType = "KO"

PCAresults_sortedSamples$sampleName = rownames(PCAresults_sortedSamples)
rownames(PCAresults_sortedSamples) = NULL
PCAresults_sortedSamples$sampleType = "sort"

#Add new samples to extracted data
t03_XfromPCA = bind_rows(t03_XfromPCA, PCAresults_sortedSamples)

#Combine metaData with PCA results
t03_XfromPCA = left_join(t03_XfromPCA, metadata, by = c("sampleName" = "geneName"))

#Fill out a few of the empty fields in the sorted samples which are needed for the plot
t03_XfromPCA$medianlFC_DP[is.na(t03_XfromPCA$medianlFC_DP)] = 2
t03_XfromPCA$medianlFC_VemR[is.na(t03_XfromPCA$medianlFC_VemR)] = 2

#Fill out the tier for sorted samples
t03_XfromPCA$DPtier[is.na(t03_XfromPCA$DPtier)] = 1
t03_XfromPCA$VemRtier[is.na(t03_XfromPCA$VemRtier)] = 1

#Fill out the effect for sorted samples
t03_XfromPCA$EffectOnScreen_DP = as.character(t03_XfromPCA$EffectOnScreen_DP)
t03_XfromPCA$EffectOnScreen_VemR = as.character(t03_XfromPCA$EffectOnScreen_VemR)
t03_XfromPCA$EffectOnScreen_DP[is.na(t03_XfromPCA$EffectOnScreen_DP)] = "sortedSample"
t03_XfromPCA$EffectOnScreen_VemR[is.na(t03_XfromPCA$EffectOnScreen_VemR)] = "sortedSample"

#Summary on PCA results
t04 = summary(PCAresults)

#save info regarding proportion of variance accounted for by each PC
t04_importance = data.frame(t04$importance)

#Change names for clarity
t03_XfromPCA$sampleName = gsub("DP", "NGFRhi/EGFRhi_Ed", t03_XfromPCA$sampleName)
t03_XfromPCA$sampleName = gsub("NoDrug_EGFR_C", "EGFRhi_Ben", t03_XfromPCA$sampleName)
t03_XfromPCA$sampleName = gsub("NoDrug_EGFRhi_B", "EGFRhi_Syd", t03_XfromPCA$sampleName)
t03_XfromPCA$sampleName = gsub("NoDrug_NGFR_C", "NGFRhi_Ben", t03_XfromPCA$sampleName)
t03_XfromPCA$sampleName = gsub("Resistant_mix_A", "ResistanceCells_IsolatedColonies_Syd", t03_XfromPCA$sampleName)
t03_XfromPCA$sampleName = gsub("Resistant_mix_B", "ResistanceCells_SortExperiment_Syd", t03_XfromPCA$sampleName)

#Obtain loadings
loadings = data.frame(PCAresults$rotation)
loadings$gene_ID = rownames(loadings)

#Plot PCA results: NGFRhi EGFRhi screen
tpmPCAplot_test <- ggplot(t03_XfromPCA, aes(PC1, PC2)) +
  geom_point(aes(color = EffectOnScreen_DP, size = (1/as.numeric(DPtier)))) +
  geom_text_repel(data = t03_XfromPCA, aes(PC1, PC2, label = sampleName, color = EffectOnScreen_DP)) +
  scale_color_manual(values = c("none" = "#CCCCCC", "up" = "#009933", "down" = "black", "sortedSample" = "#99CCFF")) +
  theme_classic()
plot(tpmPCAplot_test)
ggsave("04_plots/190321_KO-RNAseq/PCA_ESmatrix_GOBPset_PC1andPC2_byNGFRhiEGFRHiEffect_simplified.PDF", width = 10, height = 7)

#Plot PCA results: VemR screen
tpmPCAplot_test <- ggplot(t03_XfromPCA, aes(PC1, PC2)) +
  geom_point(aes(color = EffectOnScreen_VemR, size = (1/as.numeric(VemRtier)))) +
  geom_text_repel(data = t03_XfromPCA, aes(PC1, PC2, label = sampleName, color = EffectOnScreen_VemR)) +
  scale_color_manual(values = c("none" = "#CCCCCC", "up" = "#FF0000", "down" = "black", "sortedSample" = "#99CCFF")) +
  theme_classic()
plot(tpmPCAplot_test)
ggsave("04_plots/190321_KO-RNAseq/PCA_ESmatrix_GOBPset_PC1andPC2_byVemREffect_simplified.PDF", width = 10, height = 7)
