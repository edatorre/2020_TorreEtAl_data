#This plots the messages per person in slack

#load metadata
metadata = read.table("02_metadata/190824_inVivoAnalysis/190824_mouse_ID_map.txt", header = T)
metadata_measurementGroup = read.table("02_metadata/190824_inVivoAnalysis/190824_measurementGroup.txt", header = T)

#Read in and combine raw data
files <- list.files(path = "01_rawData/190824_inVivoAnalysis/tabs/", pattern=".txt")
inputData <- NULL

for (f in files) {
  dat = read.table(paste("01_rawData/190824_inVivoAnalysis/tabs/", f, sep = ""), sep = "\t", fill = TRUE, header = T)
  dat = dat %>%
    dplyr::select("Mouse_ID", "Volume_.mm3.", "Time_on_treatment_.days.")
  inputData <- rbind(inputData, dat)
}

#Clean up table: create mouse ID that matches metadata
inputData = inputData %>%
  mutate(mouse_ID = paste("mouse", "_", inputData$Mouse_ID, sep = ""))

inputData = inputData %>%
  filter(mouse_ID %in% metadata$mouse_ID) %>%
  dplyr::select(-Mouse_ID)

#Remove rows with no measurements
inputData = inputData[is.na(inputData$Time_on_treatment_.days.) == F, ]

#Remove RIPs
inputData = inputData %>%
  filter(Volume_.mm3. != "RIP") %>%
  filter(Time_on_treatment_.days. != "RIP")

#Convert character columns into numbers
inputData$Volume_.mm3. = as.numeric(inputData$Volume_.mm3.)
inputData$Time_on_treatment_.days. = as.numeric(inputData$Time_on_treatment_.days.)

#Remove new NAs from table
inputData = inputData[is.na(inputData$Time_on_treatment_.days.) == F, ]

#Add metadata to file
inputData = left_join(inputData, metadata, by = "mouse_ID") %>%
  left_join(., metadata_measurementGroup, by = "Time_on_treatment_.days.")

#Modify column names
colnames(inputData) = c("tumorVolume", "treatmentDuration", "mouseID", "KO", "treatmentArm", "measurementGroup")

#save table
write.table(inputData, file = "03_extractedData/190824_inVivoAnalysis/reformatedData.txt", quote = F, col.names = T, row.names = F, sep = "\t")


