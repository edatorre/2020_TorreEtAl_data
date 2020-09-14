#This script excludes a given mouse from the data.


excludeMouse = function(mouseToExclude) {
  
  allData = read.table("03_extractedData/190824_inVivoAnalysis/reformatedData.txt", header = T, sep = "\t")
  
  #Exclude mouse
  allData_filtered = allData %>%
    filter(mouseID != mouseToExclude)
  
  #save data
  write.table(allData_filtered, file = "03_extractedData/190824_inVivoAnalysis/reformatedData.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}


