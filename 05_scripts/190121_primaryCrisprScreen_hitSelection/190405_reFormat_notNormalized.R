# #The goal of this script is to reformat the files and merge them so that I can use them in the hit selection script. 
# 
# #Load libraries that will be used in the scripts.
# library(tidyverse)
# 
# #clear workspace
# rm(list=ls())

#If you want to hardcode the path uncomment the line below:
pathToPaperFolder = "01_rawData/190121_primaryCrisprScreen_hitSelection/rawReads/"
setwd(pathToPaperFolder)

#Read input files (reads for each of the three sgRNA libraries)
epigenetic_input = read.table("EpigeneticLibrary_Screen4/rawCounts/epigeneticScreen_combinedReads.txt", header = TRUE)
kinase_input = read.table("KinaseLibrary_Screen1/rawCounts/KinaseScreen_combinedReads.txt", header = TRUE)
TF_input = read.table("TFLibrary_Screen1/rawCounts/TFScreen_combinedReads.txt", header = TRUE)

#Read in the list of tagrets.
targets =  read.table("../../../02_metadata/ScreenTargets_metadata.txt", header = TRUE)

#Add a label indicating the library of origin to the samples.
epigenetic_input = epigenetic_input %>%
  mutate(library = "epigenetic")
kinase_input = kinase_input %>%
  mutate(library = "kinase")
TF_input = TF_input %>%
  mutate(library = "TF")

#Merge all three libraires into a single file
mergedScreen = bind_rows(epigenetic_input, kinase_input, TF_input)

#now merge the metadata file (containing information about the targets) to the raw reads file.
mergedScreen = left_join(targets, mergedScreen, by = c("sequence", "library", "gRNA"))


#save file
write.table(mergedScreen, file = "CRISPRscreens_combined.txt", sep = "\t", quote = FALSE, row.names = FALSE)








