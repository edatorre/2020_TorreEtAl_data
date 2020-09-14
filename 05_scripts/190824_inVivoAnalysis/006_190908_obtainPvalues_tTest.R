
# #For testing
# currentKO = "DOT1L"
# currentTreatmentArm = "treated"
# directionality = "greater"

runTtest = function(currentKO, currentTreatmentArm, directionality) {
  
  #Load table with volume changes
  tableOfChanges = read.table("03_extractedData/190824_inVivoAnalysis/changeInTumorVolume.txt", sep = "\t", header = TRUE)
  
  #Limit data to tretment arm
  tableOfChanges = tableOfChanges %>%
    filter(treatmentArm == currentTreatmentArm)
  
  #Make list of groups
  numberOfGroups = c("group1", "group2", "group3", "group4", "group5", "group6", "group7", "group8", "group9", "group10", "group11", "group12")
  
  #Create table to store data
  tableOfpValues = data.frame(matrix(0,length(numberOfGroups),1))
  colnames(tableOfpValues) = currentKO
  rownames(tableOfpValues) = numberOfGroups
  
  #Create index to the output table  
  group = 0
  
  for(currentGroup in numberOfGroups) {
    
    #These are my table indexes
    group = group + 1
    
    #separate controls
    control_tableOfChanges = tableOfChanges %>%
      filter(KO == "Neg") %>%
      filter(measurementGroup == currentGroup)
    
    #separate test samples
    ko_tableOfChanges = tableOfChanges %>%
      filter(KO == currentKO) %>%
      filter(measurementGroup == currentGroup)
    
    if(nrow(control_tableOfChanges) >= 3 & nrow(ko_tableOfChanges) >= 3) {
      
      #Carry out t test
      ttest = t.test(ko_tableOfChanges$logFoldChangeInTumorVolume, control_tableOfChanges$logFoldChangeInTumorVolume, alternative = directionality)
      
      #Uses the idexes above to store the p values
      tableOfpValues[group, 1] = ttest$p.value
      
    }
    
  }
  
  #save output table
  outputFileName = paste("03_extractedData/190824_inVivoAnalysis/ttests_", currentKO, "_", currentTreatmentArm,"_",directionality, ".txt", sep = "")
  write.table(tableOfpValues, file = outputFileName, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
}

