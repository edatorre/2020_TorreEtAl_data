# #Load libraries that will be used in the scripts.
# library(tidyverse)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(biomaRt)
# 
# #clear workspace
# rm(list=ls())

#Screen tagrets
screenTargets = read.table("02_metadata/ScreenTargets_metadata.txt", header = T)

epiTargets = screenTargets %>%
  filter(library == "epigenetic") %>%
  dplyr::select(gene) %>%
  distinct()

kinaseTargets = screenTargets %>%
  filter(library == "kinase") %>%
  dplyr::select(gene) %>%
  distinct()

tfTargets = screenTargets %>%
  filter(library == "TF") %>%
  dplyr::select(gene) %>%
  distinct()

epiTargets_converted = bitr(epiTargets$gene, fromType="SYMBOL", toType=("ENTREZID"), OrgDb="org.Hs.eg.db")
kinaseTargets_converted = bitr(kinaseTargets$gene, fromType="SYMBOL", toType=("ENTREZID"), OrgDb="org.Hs.eg.db")
tfTargets_converted = bitr(tfTargets$gene, fromType="SYMBOL", toType=("ENTREZID"), OrgDb="org.Hs.eg.db")

#compute enrichment analysis - upregulated genes
ego1_epi <- enrichGO(gene          = epiTargets_converted$ENTREZID,
                 #universe      = backgorundTranscriptome_converted$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
#head(ego1_epi)
ego1.df_epi = data.frame(ego1_epi)
rownames(ego1.df_epi) = NULL

ego1_kinase <- enrichGO(gene          = kinaseTargets_converted$ENTREZID,
                     #universe      = backgorundTranscriptome_converted$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
#head(ego1_kinase)
ego1.df_kinase = data.frame(ego1_kinase)
rownames(ego1.df_kinase) = NULL

ego1_tf <- enrichGO(gene          = tfTargets_converted$ENTREZID,
                        #universe      = backgorundTranscriptome_converted$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "fdr",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)
#head(ego1_tf)
ego1.df_tf = data.frame(ego1_tf)
rownames(ego1.df_tf) = NULL

write.table(ego1.df_epi, "03_extractedData/190121_primaryCrisprScreen_hitSelection/GOenrichment_epiLib.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(ego1.df_kinase, "03_extractedData/190121_primaryCrisprScreen_hitSelection/GOenrichment_kinaseLib.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
write.table(ego1.df_tf, "03_extractedData/190121_primaryCrisprScreen_hitSelection/GOenrichment_TFLib.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#Manually select GO terms to organize hits into biologically interesting sets.  













