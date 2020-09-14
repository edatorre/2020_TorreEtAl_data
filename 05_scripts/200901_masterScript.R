#This is the master script for the paper

#Clean workspace
rm(list=ls())

#NOW SET THE PATH TO THE DATA FOLDER. BELOW IS AN EXAMPLE
pathwayToData = "/Users/eduardotorre/Dropbox (RajLab)/Eduardo_shared/et_paperCRISPR/200901_data/"

#NOW JUST RUN THE REST OF THE SCRIPT...

#### Package installation ####
#This section checks if you have the packages required. If you dont, then it downlaod them. 
list.of.packages.base <- c("BiocManager","tidyverse", "stringr", "pheatmap", "RColorBrewer",
                        "ggrepel", "eulerr")
new.packages <- list.of.packages.base[!(list.of.packages.base %in% rownames(installed.packages()))]
if(length(new.packages)) install.packages(new.packages)

list.of.packages.bioconductor <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db",
                           "biomaRt", "KEGGREST", "qusage")
new.packages.bioconductor <- list.of.packages.bioconductor[!(list.of.packages.bioconductor %in% rownames(installed.packages()))]
if(length(new.packages.bioconductor)) BiocManager::install(new.packages.bioconductor)
######

#### Load libraries ####
library(tidyverse)
library(stringr)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(qusage)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(KEGGREST)
library(ggrepel)
library(eulerr)

##############################
###Primary CRISPR screen### 
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190405_combineData.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190405_reFormat_notNormalized.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190405_reFormat_normalized.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190405_obtainGOtermsOfTargets.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190405_AssignGenesToSelectedGOs_epigeneticTargets.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190405_AssignGenesToSelectedGOs_kinaseTargets.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190405_AssignGenesToSelectedGOs_tfTargets.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190429_primaryCRISPRscreen_identifyHits_function.R")
selectHits_primaryScreen("epigenetic", "DP")
selectHits_primaryScreen("kinase", "DP")
selectHits_primaryScreen("TF", "DP")
selectHits_primaryScreen("epigenetic", "VemR")
selectHits_primaryScreen("kinase", "VemR")
selectHits_primaryScreen("TF", "VemR")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190429_plotDropoutOfControls.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190402_createMetadataFile_newFormat.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/200831_calculateMedianLogFC_TargetedScreenGuides.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190121_primaryCrisprScreen_hitSelection/190416_plotOverlapBetweenPrimaryScreens.R")

##############################
### NGFR IF ###
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/180925_NGFR-IF/200831_analyzeIF.R")

##############################
### Colony growth validation of hits ###
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/181101_ColonyGrowth36_validationOfKOs_reanalyzed/190405_extractData_Plate1.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/181101_ColonyGrowth36_validationOfKOs_reanalyzed/190405_extractData_Plate2.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/181101_ColonyGrowth36_validationOfKOs_reanalyzed/190405_extractData_Plate3.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/181101_ColonyGrowth36_validationOfKOs_reanalyzed/190405_extractData_Plate4.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/181101_ColonyGrowth36_validationOfKOs_reanalyzed/190405_extractData_Plate5.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/181101_ColonyGrowth36_validationOfKOs_reanalyzed/190405_extractData_Plate6.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/181101_ColonyGrowth36_validationOfKOs_reanalyzed/190405_plotData_bothScreens.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/190430_extractData_plate1.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/190430_extractData_plate2.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/190430_extractData_plate3.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/180611_ColonyGrowth22_WM989-5a3_ScreenFollowUp_18KOs_reanalysis/190430_plotData.R")

##############################
### Validation Tables per tier ####  
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190406_validationTables/190408_createValidationTables_IFandCG.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190406_validationTables/190408_plotSlidingScale_DP.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190406_validationTables/190407_plotSlidingScale_VemR.R")

##############################
###RNA SEQUENCING OF SORTED MELANOMA CELLS#####
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190408_RNAseq_otherSamples/NGFRhiEGFRhi/190408_normalizeToRPM.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190408_RNAseq_otherSamples/NGFRhiEGFRhi/190408_runDEseq2.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190408_RNAseq_otherSamples/NGFRhiEGFRhi/190408_script1_makeGCTAndPhenotypeFiles.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190408_RNAseq_otherSamples/NGFRhiEGFRhi/190408_script2_RunGSEA.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190408_RNAseq_otherSamples/NGFRhiEGFRhi/190408_script3_createESfilesForEachDatabase.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190408_RNAseq_otherSamples/otherSamples/190408_meltCounts.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190408_RNAseq_otherSamples/otherSamples/190408_normalizeCounts.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190408_RNAseq_otherSamples/otherSamples/190408_script1_makeGCTAndPhenotypeFiles.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190408_RNAseq_otherSamples/otherSamples/190408_script2_RunGSEA.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190408_RNAseq_otherSamples/otherSamples/190408_script3_createESfilesForEachDatabase.R")

##############################
####KO RNA sequencing####
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190408_mergeRawReads.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190408_convertToRPM.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190408_plotNumCountsPerSample.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190408_runDESeq_runningCommand.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190505_clusterBasedOnDiffExpTranscriptome_1KOorMore.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190408_GSEA_script1_makeGCTAndPhenotypeFiles.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190408_GSEA_script2_RunGSEA_allKO_onLongGeneSets.R") #THIS TAKES A LONG TIME 

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190408_GSEA_script3_createESfilesForEachDatabase_onLongGeneSets.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190505_GSEA_script4_generateClustersAndHeatmaps.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
selectedGeneSetCluster = 5 #this determine which clusters to work on
selectedKOcluster = 1 #this determine which clusters to work on
source("05_scripts/190321_KO-RNAseq/190408_obtainTheLeadingEdges.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
selectedGeneSetCluster = 6 #this determine which clusters to work on
selectedKOcluster = 1 #this determine which clusters to work on
source("05_scripts/190321_KO-RNAseq/190408_obtainTheLeadingEdges.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
selectedGeneSetCluster = 8 #this determine which clusters to work on
selectedKOcluster = 1 #this determine which clusters to work on
source("05_scripts/190321_KO-RNAseq/190408_obtainTheLeadingEdges.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
selectedGeneSetCluster = 4 #this determine which clusters to work on
selectedKOcluster = 2 #this determine which clusters to work on
source("05_scripts/190321_KO-RNAseq/190408_obtainTheLeadingEdges.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190408_runClusterProfilerOnLeadingEdges.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190505_PCAbasedOnESmatrix_GOBPset_wSortedSamples_allKO_simplified.R")
source("05_scripts/190321_KO-RNAseq/190505_PCAbasedOnESmatrix_GOBPset_wSortedSamples_allKO.R")

##This is to obtain the differentiation plots
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190321_KO-RNAseq/190526_plotDifferentiationStatusVsRessitance.R")

################################
#### Targeted CRISPR screen ####
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190122_targetedCrisprScreen_hitSelection/190408_mergedRawReads.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190122_targetedCrisprScreen_hitSelection/190408_normalizeToRpm.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190122_targetedCrisprScreen_hitSelection/190409_targetedCRISPRscreen_identifyHits_function.R")
selectHits_targetedScreen("WM989.bulk", "DP")
selectHits_targetedScreen("X451Lu", "DP")
selectHits_targetedScreen("WM989.bulk", "VemR")
selectHits_targetedScreen("X451Lu", "VemR")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190122_targetedCrisprScreen_hitSelection/200831_plotScreenResults_DP.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190122_targetedCrisprScreen_hitSelection/200831_plotScreenResults_VemR.R")

################################
#Plot Immunofluorescence over colony growth ***This one has to run AFTER the KO-RNAseq scripts are run
workspace = ls() 
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190218_NGFR-IFvsResistance/200903_plotIFvsResistance copy.R")

################################
# Plot EPZ5676 timing data
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/180817_ColonyGrowth31_EPZ5676_Timing/200831_analyzeTimingData.R")

################################
#Plot cell growth effects upon EPZ5676 treatment
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/181106_cellDoublingOnDrugs_EPZ5676/190508_calculateCellDoublings_EPZ5676.R")

################################
#Obtain numbers for manuscript
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190422_obtainNumbersForManuscript/190422_obtainNumbers.R")

################################
#Create master table
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190430_makeMasterTable/190430_makeMasterTable.R")

################################
#In vivo analysis
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190824_inVivoAnalysis/001_190825_formatData.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190824_inVivoAnalysis/001-1_excludeMouse.R")
excludeMouse("mouse_39") #From BRD2 treated arm. Large deviations because two tumors merged.
excludeMouse("mouse_46") #The change in tumor volume between day 14 and day 21 is suggestive of a technical error. 

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190824_inVivoAnalysis/003_190825_obtainChangeInTumorVolume.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190824_inVivoAnalysis/004_190825_plotSummaryStats_truncatedPerArm.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190824_inVivoAnalysis/006_190908_obtainPvalues_tTest.R")

runTtest("DOT1L", "treated", "greater")
runTtest("DOT1L", "untreated", "greater")
runTtest("DOT1L", "untreated", "less")

runTtest("LATS2", "treated", "greater")
runTtest("LATS2", "untreated", "greater")
runTtest("LATS2", "untreated", "less")

runTtest("BRD2", "treated", "less")
runTtest("BRD2", "untreated", "less")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/190824_inVivoAnalysis/007_200901_nMiceAndPvals.R")

################################
#Comparison of resistance profile betwen WM989 and WM989-cas9 
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/200302_PLXdosing_limitedSet_WM989comp/200831_compareCellLines_colonyGrowth.R")

################################
#Comparison of KOs at different concentrations of PLX4032 
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/200302_PLXdosing_limitedSet_responseCurve/200831_PLXresponseCurve_colonyGrowth.R")

################################
#Transcriptome profiling of subpopulations of EPZ5676-treated and DMSO-treated melanoma cells
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/200309_EPZtreatedCells_RNAseqData/200330_renameInitialFile.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/200309_EPZtreatedCells_RNAseqData/200330_convertToRPM.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/200309_EPZtreatedCells_RNAseqData/200330_plotNumCountsPerSample.R")

workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/200309_EPZtreatedCells_RNAseqData/200330_PCAonGeneExpession.R")

################################
#Colony growth on subpopulations of melanoma cells after EPZ5676 treatment
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/200707_EPZtreatment_CJET/200831_plotColonyGrowth_subpopulations.R")

################################
#NegativeSelectionAssay for Cas9 cell lines
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/170215_WM989cas9_NegativeSelectionAssay/190311_negativeSelectionAssay_plotData.R")

################################
#Colony growth for subpopulations of WM989-cas9 wihtout any treatment
workspace = ls()
workspace = workspace[workspace != "pathwayToData"]
rm(list = workspace)
setwd(pathwayToData)
source("05_scripts/170731_Cas9_ColonyGrowthAssays/190311_plotData.R")

