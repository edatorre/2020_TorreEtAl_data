# #This script compiles the reads for each of the samples and each of the libraries. 
# #I originally had three sgRNA libraries (TF, kinase, and epigenetic libraries). The code below is split into three section. 
# #In each section I work with a different library. 


#If you want to hardcode the path uncomment the line below:
pathToPaperFolder = "01_rawData/190121_primaryCrisprScreen_hitSelection/rawReads/"
setwd(pathToPaperFolder)

#Combine epigenetic library reads.
#Load datasets and modify column names

  #Seq1
  epi_baseline_seq1 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq1/8a_Baseline1_Epi.txt")
  colnames(epi_baseline_seq1) = c("gRNA", "sequence", "baseline_1") 
  epi_ungated_seq1 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq1/7a_Baseline2_Epi (2).txt")
  colnames(epi_ungated_seq1) = c("gRNA", "sequence", "ungated_1")
  epi_DP_seq1 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq1/7d_DP_Epi (1).txt")
  colnames(epi_DP_seq1) = c("gRNA", "sequence", "DP_1")
  epi_VemR_seq1 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq1/Epi Lib #4 d64 Vem. Resistant 1000x.txt")
  colnames(epi_VemR_seq1) = c("gRNA", "sequence", "VemR_1")
  epi_NGFR_seq1 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq1/7b_NGFR_Epi (1).txt")
  colnames(epi_NGFR_seq1) = c("gRNA", "sequence", "NGFR_1")
  epi_EGFR_seq1 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq1/7c_EGFR_Epi (1).txt")
  colnames(epi_EGFR_seq1) = c("gRNA", "sequence", "EGFR_1")
  
  #Seq2
  epi_baseline_seq2 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq2/Epi Baseline 1_3.txt")
  colnames(epi_baseline_seq2) = c("gRNA", "sequence", "baseline_2")
  epi_ungated_seq2 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq2/Epi Baseline 2_3.txt")
  colnames(epi_ungated_seq2) = c("gRNA", "sequence", "ungated_2")
  epi_DP_seq2 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq2/Epi DP_3.txt")
  colnames(epi_DP_seq2) = c("gRNA", "sequence", "DP_2")
  epi_NGFR_seq2 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq2/Epi NGFR+_3.txt")
  colnames(epi_NGFR_seq2) = c("gRNA", "sequence", "NGFR_2")
  epi_EGFR_seq2 = read.table("EpigeneticLibrary_Screen4/rawCounts/seq2/Epi EGFR+_3.txt")
  colnames(epi_EGFR_seq2) = c("gRNA", "sequence", "EGFR_2")
  
#Combine all samples into a single file
  t01_allData = left_join(epi_baseline_seq1, epi_baseline_seq2, by = c("gRNA", "sequence"))%>%
    left_join(., epi_ungated_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., epi_ungated_seq2, by = c("gRNA", "sequence")) %>%
    left_join(., epi_DP_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., epi_DP_seq2, by = c("gRNA", "sequence")) %>%
    left_join(., epi_VemR_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., epi_NGFR_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., epi_NGFR_seq2, by = c("gRNA", "sequence")) %>%
    left_join(., epi_EGFR_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., epi_EGFR_seq2, by = c("gRNA", "sequence"))
  
#Create file where I combine (by addition) the reads from multiple sequencing runs  
  t02_jointData = t01_allData %>%
    dplyr::select(gRNA, sequence)
  t02_jointData = t02_jointData %>%
    mutate(baseline = t01_allData$baseline_1 + t01_allData$baseline_2) %>%
    mutate(ungated = t01_allData$ungated_1 + t01_allData$ungated_2) %>%
    mutate(DP = t01_allData$DP_1 + t01_allData$DP_2) %>%
    mutate(VemR = t01_allData$VemR_1) %>%
    mutate(EGFR = t01_allData$EGFR_1 + t01_allData$EGFR_2) %>%
    mutate(NGFR = t01_allData$NGFR_1 + t01_allData$NGFR_2)
  
#save table of raw values
write.table(t02_jointData, file = "EpigeneticLibrary_Screen4/rawCounts/epigeneticScreen_combinedReads.txt", quote = FALSE, sep = "\t", row.names = FALSE)

#clear workspace
rm(t01_allData)

#Repeat the process above but now for the kinase library
#Combine kinase library reads.

#Load datasets
  #Seq3 - baseline only
  kinase_baseline_seq1 = read.table("KinaseLibrary_Screen1/rawCounts/Seq3/Eduardo Kinase Baseline d6.txt")
  colnames(kinase_baseline_seq1) = c("gRNA", "sequence", "baseline_1")
  #Seq1
  kinase_ungated_seq1 = read.table("KinaseLibrary_Screen1/rawCounts/Seq1/Kinase Lib #1 d30 Baseline.txt")
  colnames(kinase_ungated_seq1) = c("gRNA", "sequence", "ungated_1")
  kinase_DP_seq1 = read.table("KinaseLibrary_Screen1/rawCounts/Seq1/Kinase Lib #1 d30 DP.txt")
  colnames(kinase_DP_seq1) = c("gRNA", "sequence", "DP_1")
  kinase_VemR_seq1_1 = read.table("KinaseLibrary_Screen1/rawCounts/Seq1/Kinase Screen B Vem. R1.txt")
  colnames(kinase_VemR_seq1_1) = c("gRNA", "sequence", "VemR_1_1")
  kinase_VemR_seq1_2 = read.table("KinaseLibrary_Screen1/rawCounts/Seq1/Kinase Screen B Vem. R2.txt")
  colnames(kinase_VemR_seq1_2) = c("gRNA", "sequence", "VemR_1_2")
  kinase_VemR_seq1_3 = read.table("KinaseLibrary_Screen1/rawCounts/Seq1/Kinase Screen B Vem. R3.txt")
  colnames(kinase_VemR_seq1_3) = c("gRNA", "sequence", "VemR_1_3")
  kinase_VemR_seq1_4 = read.table("./KinaseLibrary_Screen1/rawCounts/Seq1/Kinase Screen B Vem. R4.txt")
  colnames(kinase_VemR_seq1_4) = c("gRNA", "sequence", "VemR_1_4")
  kinase_NGFR_seq1 = read.table("KinaseLibrary_Screen1/rawCounts/Seq1/Kinase Lib #1 d30 NGFR+.txt")
  colnames(kinase_NGFR_seq1) = c("gRNA", "sequence", "NGFR_1")
  kinase_EGFR_seq1 = read.table("KinaseLibrary_Screen1/rawCounts/Seq1/Kinase Lib #1 d30 EGFR+.txt")
  colnames(kinase_EGFR_seq1) = c("gRNA", "sequence", "EGFR_1")
  #Seq2
  kinase_ungated_seq2 = read.table("KinaseLibrary_Screen1/rawCounts/Seq2/Kinase Lib #1 d30 Baseline 1000x.txt")
  colnames(kinase_ungated_seq2) = c("gRNA", "sequence", "ungated_2")
  kinase_DP_seq2 = read.table("KinaseLibrary_Screen1/rawCounts/Seq2/Kinase Lib #1 d30 DP 100x.txt")
  colnames(kinase_DP_seq2) = c("gRNA", "sequence", "DP_2")
  kinase_NGFR_seq2 = read.table("KinaseLibrary_Screen1/rawCounts/Seq2/Kinase Lib #1 d30 NGFR+ 1000x.txt")
  colnames(kinase_NGFR_seq2) = c("gRNA", "sequence", "NGFR_2")
  kinase_EGFR_seq2 = read.table("KinaseLibrary_Screen1/rawCounts/Seq2/Kinase Lib #1 d30 EGFR+ 100x.txt")
  colnames(kinase_EGFR_seq2) = c("gRNA", "sequence", "EGFR_2")
  
  
  #Combine all samples into a single file
  t01_allData = left_join(kinase_baseline_seq1, kinase_ungated_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_ungated_seq2, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_DP_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_DP_seq2, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_VemR_seq1_1, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_VemR_seq1_2, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_VemR_seq1_3, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_VemR_seq1_4, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_NGFR_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_NGFR_seq2, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_EGFR_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., kinase_EGFR_seq2, by = c("gRNA", "sequence"))
  
  t02_jointData = t01_allData %>%
    dplyr::select(gRNA, sequence) 
  
  t02_jointData = t02_jointData %>%
    mutate(baseline = t01_allData$baseline_1) %>%
    mutate(ungated = t01_allData$ungated_1 + t01_allData$ungated_2) %>%
    mutate(DP = t01_allData$DP_1 + t01_allData$DP_2) %>%
    mutate(VemR = t01_allData$VemR_1_1 + t01_allData$VemR_1_2 + t01_allData$VemR_1_3  + t01_allData$VemR_1_4) %>%
    mutate(NGFR = t01_allData$NGFR_1 + t01_allData$NGFR_2) %>%
    mutate(EGFR = t01_allData$EGFR_1 + t01_allData$EGFR_2)
  
  #save table of raw values
  write.table(t02_jointData, file = "KinaseLibrary_Screen1/rawCounts/KinaseScreen_combinedReads.txt", quote = FALSE, sep = "\t", row.names = FALSE)
  
  #clear workspace
  rm(t01_allData)

#Repeat the process above but now for the TF library 
#Combine TF library reads.
  
  #Load datasets
  TF_baseline_seq1 = read.table("./TFLibrary_Screen1/rawCounts/Eduardo TF Baseline d6.txt")
  colnames(TF_baseline_seq1) = c("gRNA", "sequence", "baseline")
  TF_ungated_seq1 = read.table("./TFLibrary_Screen1/rawCounts/TF Lib #1 d30 Baseline 1000x.txt")
  colnames(TF_ungated_seq1) = c("gRNA", "sequence", "ungated")
  TF_DP_seq1 = read.table("./TFLibrary_Screen1/rawCounts/TF #1 d30 DP 75x.txt")
  colnames(TF_DP_seq1) = c("gRNA", "sequence", "DP")
  TF_VemR_seq1 = read.table("./TFLibrary_Screen1/rawCounts/TF Screen B Vem. R.txt")
  colnames(TF_VemR_seq1) = c("gRNA", "sequence", "VemR")
  TF_NGFR_seq1 = read.table("./TFLibrary_Screen1/rawCounts/TF #1 d30 NGFR+ 500x.txt")
  colnames(TF_NGFR_seq1) = c("gRNA", "sequence", "NGFR")
  TF_EGFR_seq1 = read.table("./TFLibrary_Screen1/rawCounts/TF #1 d30 EGFR+ 170x.txt")
  colnames(TF_EGFR_seq1) = c("gRNA", "sequence", "EGFR")
  
  tmp = TF_baseline_seq1 %>%
    dplyr::select(gRNA, sequence)
  
  t01_allData = left_join(TF_baseline_seq1, TF_ungated_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., TF_DP_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., TF_VemR_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., TF_NGFR_seq1, by = c("gRNA", "sequence")) %>%
    left_join(., TF_EGFR_seq1, by = c("gRNA", "sequence"))
  
  t02_jointData = left_join(tmp, t01_allData, by = c("gRNA", "sequence"))
  
  #save table of raw values
  write.table(t02_jointData, file = "TFLibrary_Screen1/rawCounts/TFScreen_combinedReads.txt", quote = FALSE, sep = "\t", row.names = FALSE)

  
  
    