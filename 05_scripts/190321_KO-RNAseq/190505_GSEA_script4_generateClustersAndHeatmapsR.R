# #This script is to cluster KOs by enrichment score of gene sets. It also generates a heatmap of those results. 

#Detemrine pathway to ES.combined files
path.to.ES.files = "03_extractedData/190321_KO-RNAseq/GSEA/ESfiles/"

#Select controls to eliminate from analysis
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
negativeControls = c("neg", "WT")

#Load metadata
metadata = read.table("02_metadata/190402_metadataOfHits.txt", header = T)
metadata = metadata %>%
  filter(!geneName %in% positiveControls) %>%
  filter(!geneName %in% negativeControls)
metadata$geneName = gsub("-", "\\.", metadata$geneName)

  #Make list of ES files in directory
  ES.files = list.files(path = path.to.ES.files, full.names = FALSE)

  #Select the gene set here!!!! Just change the number
  combinedFile = ES.files[1]
 
  print(paste("working on ", combinedFile, sep = ""))
  flush.console()
  
  #Read in enrichment score matrix
  enrichmentScoreMatrix = read.table(paste(path.to.ES.files, combinedFile, sep = ""), header = T)
  enrichmentScoreMatrix$GS = as.character(enrichmentScoreMatrix$GS)
  
  #Read in enrichment score of new samples: NGFRhi EGFRhi and VemR
  enrichmentScoreMatrix_DP = read.table("03_extractedData/190408_RNAseq_otherSamples/NGFRandEGFRhi/GSEA/ESfiles/c5.bp.v6.2_combinedESfile.txt", header = T)
  enrichmentScoreMatrix_DP$GS = as.character(enrichmentScoreMatrix_DP$GS)
  enrichmentScoreMatrix_otherSamples = read.table("03_extractedData/190408_RNAseq_otherSamples/otherDatatsets/GSEA/ESfiles/c5.bp.v6.2_combinedESfile.txt", header = T)
  enrichmentScoreMatrix_otherSamples$GS = as.character(enrichmentScoreMatrix_otherSamples$GS)
  
  #Merge the two new files into a single one
  enrichmentScoreMatrix_newSamples = left_join(enrichmentScoreMatrix_DP, enrichmentScoreMatrix_otherSamples, by = "GS")
  
  #Change names for Arjun
  colnames(enrichmentScoreMatrix_newSamples) = c("GS", "NGFRhi/EGFRhi_Ed", "EGFRhi_Ben", "EGFRhi_Syd", "NGFRhi_Ben", "ResistanceCells_IsolatedColonies_Syd", "ResistanceCells_SortExperiment_Syd")
  
  #Keep only rows with enrichment scores > 0 for at least one KO
  enrichmentScoreMatrix = enrichmentScoreMatrix[rowSums(abs(enrichmentScoreMatrix[2:ncol(enrichmentScoreMatrix)]) != 0) >= 1, ]
  
  #Remove row with NA
  enrichmentScoreMatrix = enrichmentScoreMatrix[!is.na(enrichmentScoreMatrix$GS) ,]
  
  #Remove column with NA
  enrichmentScoreMatrix = enrichmentScoreMatrix[ , colSums(is.na(enrichmentScoreMatrix)) == 0]
  
  #Convert first column into row names
  rownames(enrichmentScoreMatrix) = enrichmentScoreMatrix$GS
  enrichmentScoreMatrix$GS = NULL
  
  #Eliminate unwated KOs from table
  enrichmentScoreMatrix = enrichmentScoreMatrix %>%
    dplyr::select(one_of(as.character(metadata$geneName)))
  
  #Now merge Data
  enrichmentScoreMatrix$GS = rownames(enrichmentScoreMatrix)
  enrichmentScoreMatrix = left_join(enrichmentScoreMatrix, enrichmentScoreMatrix_newSamples, by = "GS")
  rownames(enrichmentScoreMatrix) = enrichmentScoreMatrix$GS
  enrichmentScoreMatrix$GS = NULL
  
  #Convert table into matrix
  enrichmentScoreMatrix.matrix = as.matrix(t(enrichmentScoreMatrix))
  
  #Genearte color palete
  color4heatmap <- rev(brewer.pal(11,'RdYlBu'))
  
  color4heatmap_mid <- c('#1e215c','#2b2f82', color4heatmap,
                         '#590014','#270009')
  
  color4heatmap_dark <- c('#0b0c22','#1e215c','#2b2f82', color4heatmap,
                          '#590014','#270009','#0d0003')
  
  #create name of the plot
  databaseName = strsplit(combinedFile, "_")
  
  #create annotations table to heatmap
  KO_annotations = metadata 
  KO_annotations = KO_annotations %>%
    dplyr::select(geneName, medianlFC_DP, medianlFC_VemR, EffectOnScreen_DP, EffectOnScreen_VemR, DPtier, VemRtier) %>%
    distinct()
  KO_annotations = data.frame(KO_annotations)
  KO_annotations = KO_annotations %>%
    filter(geneName %in% rownames(enrichmentScoreMatrix.matrix)) %>%
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

  #Make heatmap
  heatmapPlot = pheatmap(
    mat = (enrichmentScoreMatrix.matrix),
    clustering_method = "ward.D2", 
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    scale = "none",
    color = color4heatmap_mid,
    annotation_row = KO_annotations,
    annotation_colors = ann_colors,
    border_color = NA,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize = 4, 
    main = "KO heatmap - GSEA on GO terms",
    filename = "04_plots/190321_KO-RNAseq/KOheatmap_GSEAonGOterms_eucledian_wardD2.pdf"
  )

  columnCuts = 10
  rowCuts = 6
  
  #Make heatmap
  heatmapPlot = pheatmap(
    mat = (enrichmentScoreMatrix.matrix),
    clustering_method = "ward.D2", 
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    scale = "none",
    color = color4heatmap_mid,
    annotation_row = KO_annotations,
    annotation_colors = ann_colors,
    cutree_rows = rowCuts,
    border_color = NA,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize = 4, 
    main = "KO heatmap - GSEA on GO terms",
    filename = "04_plots/190321_KO-RNAseq/KOheatmap_GSEAonGOterms_eucledian_wardD2_cutRows.pdf"
  )
  
  #Make heatmap
  heatmapPlot = pheatmap(
    mat = (enrichmentScoreMatrix.matrix),
    clustering_method = "ward.D2", 
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    scale = "none",
    color = color4heatmap_mid,
    annotation_row = KO_annotations,
    annotation_colors = ann_colors,
    cutree_rows = rowCuts,
    cutree_cols = columnCuts,
    border_color = NA,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize = 4, 
    main = "KO heatmap - GSEA on GO terms",
    filename = "04_plots/190321_KO-RNAseq/KOheatmap_GSEAonGOterms_eucledian_wardD2_cutRowsAndColumns.pdf"
  )
  
  #Now obtain list of KOs that cluster together
  row_dividedTree = cutree(heatmapPlot$tree_row, k = rowCuts)
  row_dividedTree = data.frame(row_dividedTree)
  row_dividedTree$KO = rownames(row_dividedTree)
  row_dividedTree = row_dividedTree[heatmapPlot$tree_row$order,]
  
  #Now obtain list of Gene sets that cluster together
  col_dividedTree = cutree(heatmapPlot$tree_col, k = columnCuts)
  col_dividedTree = data.frame(col_dividedTree)
  col_dividedTree$GS = rownames(col_dividedTree)
  col_dividedTree = col_dividedTree[heatmapPlot$tree_col$order,]

  clusterContentsFileName = paste("03_extractedData/190321_KO-RNAseq/GSEA/clusterContents/", databaseName[[1]][1], "_clusterContents_KOs.txt", sep = "")
  write.table(row_dividedTree, file = clusterContentsFileName, sep = "\t", quote = FALSE, 
              col.names = T, row.names = TRUE)
    
  clusterContentsFileName = paste("03_extractedData/190321_KO-RNAseq/GSEA/clusterContents/", databaseName[[1]][1], "_clusterContents_geneSets.txt", sep = "")
  write.table(col_dividedTree, file = clusterContentsFileName, sep = "\t", quote = FALSE, 
              col.names = T, row.names = TRUE)
  


