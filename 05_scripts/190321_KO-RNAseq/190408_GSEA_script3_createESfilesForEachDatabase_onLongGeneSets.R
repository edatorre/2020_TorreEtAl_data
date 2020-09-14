# #The goal of this script is to create a matrix of enrichment values for all my KOs and to cluster samples based on those values

#Load list of databases
listOfDatabases = list.dirs(path = "03_extractedData/190321_KO-RNAseq/GSEA/GSEAoutput/",
                            recursive = FALSE, full.names = FALSE)

#Set the path to the database gene sets
database.genesets = "02_metadata/190321_KO-RNAseq/GSEAgeneSets/"

#For each data base 
for (database in listOfDatabases) {
  
  # #for testing
  # database = listOfDatabases[1]
  
  #set the directory path to the database enrichment scores
  database.ES.path = paste("03_extractedData/190321_KO-RNAseq/GSEA/GSEAoutput/", 
                           database, "/", sep = "")
  
  #print the database 
  print(paste("working on ", database, sep = ""))
  flush.console()
  
  #Obtain gene sets
  geneSetFile = read.gmt((paste(database.genesets, database, ".symbols.gmt", sep = "")))
  geneSetFile_table = data.frame(names(geneSetFile)) 
  colnames(geneSetFile_table) = "GS"
  
  #Read the name of the KOs present
  listOfKO = list.dirs(path = database.ES.path, recursive = FALSE, full.names = FALSE)
  
  #For each KO, add the ESs to the main summary file
  for(KO in listOfKO) {
    
    #for testing
    #KO = "MITF"
    
    print(KO)
    flush.console()
    
    #Obtain file names
    path.to.ES.files = paste(database.ES.path, KO, "/", sep = "")
    files = list.files(path = path.to.ES.files, pattern = "SUMMARY.RESULTS.REPORT")
    
    #Determine how many summary files are in the folder
    numberOfFiles = length(files)
    
    #create file to store output
    GSEAfile_combined = data.frame()
    
      #Combine the values from all files for a given KO
      if (numberOfFiles > 0) {
        
          #Read in summary files
          for(j in files) {
            #for testing
            #j = files[2]
            
            #read in the file
            GSEAfile = read.table(paste(path.to.ES.files, j, sep = ""), header = T, sep = "\t")
            
            #add te contents to the combined file
            GSEAfile_combined = bind_rows(GSEAfile_combined, GSEAfile)
            
          }
        
        #Select Fields of interes only 
        GSEAfile_combined = GSEAfile_combined %>%
          dplyr::select(GS, NES)
        #Modify column name to add KO name
        colnames(GSEAfile_combined)[2] = KO
        
        #Write results to main output file
        geneSetFile_table = left_join(geneSetFile_table, GSEAfile_combined, by = "GS")
        geneSetFile_table = geneSetFile_table %>%
          distinct()
        
        #erase the GSEA file
        rm(GSEAfile_combined)
        
      }
      
    
  
  }
  
  #Convert NA's to 0 in geneSetFile_table
  geneSetFile_table[2:ncol(geneSetFile_table)][is.na(geneSetFile_table[2:ncol(geneSetFile_table)])] = 0
  
  # #Remove rows in total enrichment scores of 0
  # tmp = geneSetFile_table[rowSums(abs(geneSetFile_table[2:ncol(geneSetFile_table)]) != 0), ]
  
  #save file
  path.to.combined.ES = paste("03_extractedData/190321_KO-RNAseq/GSEA/ESfiles/", database,
                              "_combinedESfile", ".txt", sep = "")
  write.table(geneSetFile_table, file = path.to.combined.ES, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
}
