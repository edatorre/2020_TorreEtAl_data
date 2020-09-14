# #This script generates a heatmap plot of KOs using only genes Diff. Exp on NGFRhi EGFR hi cells

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

#Load ddifferential expression table
KOdata = read.table("03_extractedData/190321_KO-RNAseq/DEseq/differentialExpression_DESeq_allTargets.txt", header = T)

#convert data types
KOdata$log2FoldChange = as.numeric(as.character(KOdata$log2FoldChange))
KOdata$sampleKO = as.character(KOdata$sampleKO)
KOdata$id = as.character(KOdata$id)
KOdata$padj = as.numeric(as.character(KOdata$padj))
KOdata$baseMean = as.numeric(as.character(KOdata$baseMean))

#Remove NA's and keep unique entries
KOdata[is.na(KOdata)] = 0 
KOdata = KOdata %>%
  distinct()

#Make list of differentially expressed genes
diffExpGenes = KOdata %>%
  filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>%
  filter(padj <= 0.05) %>%
  filter(baseMean >= 30) %>%  
  group_by(id) %>%
  mutate(numberOfKO = n_distinct(sampleKO)) %>%
  ungroup(id)
  
###HERE IS WHERE I SELECT HOW MANY KOs NEED TO HAVE THE GENE DIFF. EXP###
diffExpGenes = diffExpGenes %>%
  filter(numberOfKO >= 1) %>%
  dplyr::select(id) %>%
  distinct()
###

#Limit the matrix to only genes differentially expressed
KOdata = KOdata %>%
  filter(id %in% diffExpGenes$id)

#Eliminate controls from the matrix
KOdata = KOdata %>%
  filter(!sampleKO %in% positiveControls) %>%
  filter(!sampleKO %in% negativeControls)

#Convert tall table into matrix
KOdata_matrix = KOdata %>%
  dplyr::select(id, sampleKO, log2FoldChange) %>%
  spread(sampleKO, log2FoldChange)

#Add rownames and elinimate ID field
rownames(KOdata_matrix) = KOdata_matrix$id
KOdata_matrix$id = NULL

# Remove rows and columns without all zero values  
KOdata_matrix = KOdata_matrix[, colSums(abs(KOdata_matrix) != 0) > 0]
KOdata_matrix = KOdata_matrix[rowSums(abs(KOdata_matrix[3:ncol(KOdata_matrix)])) > 0, ]

#Transpose the matrix.
KOdata_matrix_tranposed = t(KOdata_matrix)

#Create a correlation matrix
corMat <- cor(KOdata_matrix)

#Genearte color paletes
color4heatmap <- rev(brewer.pal(6,'RdYlBu'))

color4heatmap_dark <- c('#0b0c22','#1e215c','#2b2f82', color4heatmap,
                        '#590014','#270009','#0d0003')

#create annotations table to heatmap
KO_annotations = metadata 
KO_annotations = KO_annotations %>%
  dplyr::select(geneName, medianlFC_DP, medianlFC_VemR, EffectOnScreen_DP, EffectOnScreen_VemR, DPtier, VemRtier) %>%
  distinct()
KO_annotations = data.frame(KO_annotations)
KO_annotations = KO_annotations %>%
  filter(geneName %in% rownames(KOdata_matrix_tranposed)) %>%
  filter(!(geneName == "EP300" & DPtier == "tier066")) %>% #this helps me remove duplicated rows due to metadata of different 
  filter(!(geneName == "PBRM1" & DPtier == "tierNotHit")) #this helps me remove duplicated rows due to metadata of different
#filter(!(geneName == "PBRM1" & DPtier == "tier050")) #this helps me remove duplicated rows due to metadata of different
rownames(KO_annotations) = KO_annotations$geneName
KO_annotations$geneName = NULL

#Convert tiers into numbers to use later
KO_annotations$DPtier = gsub("tier075", 1, KO_annotations$DPtier)
KO_annotations$DPtier = gsub("tier066", 2, KO_annotations$DPtier)
KO_annotations$DPtier = gsub("tier050", 3, KO_annotations$DPtier)
KO_annotations$DPtier = gsub("tierNotHit", 4, KO_annotations$DPtier)

KO_annotations$VemRtier = gsub("tier075", 1, KO_annotations$VemRtier)
KO_annotations$VemRtier = gsub("tier066", 2, KO_annotations$VemRtier)
KO_annotations$VemRtier = gsub("tier050", 3, KO_annotations$VemRtier)
KO_annotations$VemRtier = gsub("tierNotHit", 4, KO_annotations$VemRtier)

#Detemine color of label
KO_annotations$labelColor = "none"

#Determine colors for labels (I select the color based on the effect on the screen where the targer was the highest tier. If there is a tie, I always select attribute the targets to the screen where it showed the largest effect)
#Select the effect of the higest tier
KO_annotations$labelColor = ifelse(KO_annotations$DPtier <= KO_annotations$VemRtier & KO_annotations$EffectOnScreen_DP == "up", "green_state", KO_annotations$labelColor)
KO_annotations$labelColor = ifelse(KO_annotations$DPtier <= KO_annotations$VemRtier & KO_annotations$EffectOnScreen_DP == "down", "gray_state", KO_annotations$labelColor)
KO_annotations$labelColor = ifelse(KO_annotations$VemRtier < KO_annotations$DPtier & KO_annotations$EffectOnScreen_VemR == "up", "red_fate", KO_annotations$labelColor)
KO_annotations$labelColor = ifelse(KO_annotations$VemRtier < KO_annotations$DPtier & KO_annotations$EffectOnScreen_VemR == "down", "gray_fate", KO_annotations$labelColor)
#When tehre is a tie in tiers, select the strongest effect. 
KO_annotations$labelColor = ifelse((KO_annotations$DPtier == KO_annotations$VemRtier) & (abs(KO_annotations$medianlFC_DP) > abs(KO_annotations$medianlFC_VemR)) &  KO_annotations$medianlFC_DP > 0, "green_state", KO_annotations$labelColor)
KO_annotations$labelColor = ifelse((KO_annotations$DPtier == KO_annotations$VemRtier) & (abs(KO_annotations$medianlFC_DP) > abs(KO_annotations$medianlFC_VemR)) &  KO_annotations$medianlFC_DP < 0, "gray_state", KO_annotations$labelColor)
KO_annotations$labelColor = ifelse((KO_annotations$DPtier == KO_annotations$VemRtier) & (abs(KO_annotations$medianlFC_DP) < abs(KO_annotations$medianlFC_VemR)) &  KO_annotations$medianlFC_VemR > 0, "red_fate", KO_annotations$labelColor)
KO_annotations$labelColor = ifelse((KO_annotations$DPtier == KO_annotations$VemRtier) & (abs(KO_annotations$medianlFC_DP) < abs(KO_annotations$medianlFC_VemR)) &  KO_annotations$medianlFC_VemR < 0, "gray_fate", KO_annotations$labelColor)

#remove extra columns
KO_annotations$medianlFC_DP = NULL
KO_annotations$medianlFC_VemR = NULL
KO_annotations$EffectOnScreen_DP = NULL
KO_annotations$EffectOnScreen_VemR = NULL

#detemrine colors for KO annotations
ann_colors = list(labelColor = c("green_state" = "green", "red_fate" = "red", "gray_state" = "gray", "gray_fate" = "gray"),
                  DPtier = c("1" = "#FFD700", "2" = "#00FFFF", "3" = "#FF00FF", "4" = "#D2691E"), 
                  VemRtier = c("1" = "#FFD700", "2" = "#00FFFF", "3" = "#FF00FF", "4" = "#D2691E")) 


#set breaks for correlation table
breaklist = seq(-3, 3, by = 0.5)

#Make heatmap
heatmapPlot = pheatmap(
  mat = as.matrix(KOdata_matrix_tranposed),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  scale = "none",
  color = color4heatmap_dark,
  annotation_row = KO_annotations,
  annotation_colors = ann_colors,
  breaks = breaklist,
  border_color = NA,
  show_colnames = FALSE,
  show_rownames = TRUE,
  filename = "04_plots/190321_KO-RNAseq/KOpheatmap_DiffExpTranscriptome_1orMoreKO.pdf",
  fontsize = 4,
  main = "KO heatmap - Diff. Exp. Transcriptome"
)

