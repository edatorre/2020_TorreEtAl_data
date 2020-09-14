### Combine reads from follow up screens ####

#Load datasets
  #Obtain sample names
    listOfFiles = list.files("01_rawData/190122_targetedCrisprScreen_hitSelection/rawCounts_modifiedNames/individualFiles/")
  #Obtain list of targets
    listOfTargets = read.table("02_metadata/TargetedScreen_targets.txt", header = TRUE)
    listOfTargets = listOfTargets %>%
      dplyr::select(gRNA)
  #Create file for merge
    rawReads_combined = listOfTargets
  #Create loop to read and merge files
    for (i in listOfFiles) {
      
      # #For testing
      # i = listOfFiles[1]
      
      #Define file name
      fileName = paste("01_rawData/190122_targetedCrisprScreen_hitSelection/rawCounts_modifiedNames/individualFiles/", i, sep = "")
     
      #read file
      inputFile = read.table(fileName, header = FALSE)
      #rename Column with reads
      colnames(inputFile)[3] = i
      #remove sequence column
      inputFile = inputFile %>%
        dplyr::select(-V2)
      #modify column names
      colnames(inputFile)[1] = "gRNA"
      
      #merge file
      rawReads_combined = left_join(rawReads_combined, inputFile, by = "gRNA")
      
    }
  

  #save combined raw reads
    write.table(rawReads_combined, file = "01_rawData/190122_targetedCrisprScreen_hitSelection/rawCounts_modifiedNames/rawReads_combined.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
