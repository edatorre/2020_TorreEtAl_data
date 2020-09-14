#the goal of this script is to calclate the differences in growth for all mice at all timepoints

#load raw data
allData = read.table("03_extractedData/190824_inVivoAnalysis/reformatedData.txt", header = T)

#Filter inittial tumor size
initialTumors = allData %>%
  filter(measurementGroup == "group1")

#Keep relevant data
initialTumors = initialTumors %>%
  dplyr::select(KO, treatmentArm, mouseID, tumorVolume)

#Rename tumor size to initial tumor size
colnames(initialTumors)[4] = "inititalTumorVolume"

#Add data to main data table
allData = left_join(allData, initialTumors, by = c("KO", "treatmentArm", "mouseID"))

#Calculate fold change from t0
allData = allData %>%
  mutate(logFoldChangeInTumorVolume = log2(tumorVolume / inititalTumorVolume))

#Save output
write.table(allData, file = "03_extractedData/190824_inVivoAnalysis/changeInTumorVolume.txt", col.names = TRUE, row.names = FALSE, quote = F, sep = "\t")
