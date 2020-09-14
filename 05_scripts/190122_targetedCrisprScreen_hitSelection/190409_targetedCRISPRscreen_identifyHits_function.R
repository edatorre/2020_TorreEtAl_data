# #This function selects hits from the screen and plots the results. 

selectHits_targetedScreen = function(CellLine, BiologicalSampleName) {

  # #For testing
  # CellLine = "X451Lu"
  # BiologicalSampleName = "DP"
  
  #Determine variables
  SelectedCellLine = CellLine
  currentSample = BiologicalSampleName
  SelectedCellLine_ungated = paste(CellLine, "_Ungated.txt", sep = "")
  SelectedCellLine_testSample = paste(CellLine, "_", currentSample ,".txt", sep = "")
  
  #Define the positive controls
  positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
  negativeControls = c("neg", "WT")
  
  #Read in metadata 
  t01_metadata = read.table("02_metadata/TargetedScreen_targets.txt", sep = "\t", header = T)
  
  #Read screen data. These are the rpm-normalized read counts for each guide and well as the unnormalized values.
  t01_allSamples_inputData = read.table("01_rawData/190122_targetedCrisprScreen_hitSelection/rawCounts_modifiedNames/normalizedReads_combined.txt", header = T, sep = "\t") #normalized
  t01_allSamples_inputData_rawCounts = read.table("01_rawData/190122_targetedCrisprScreen_hitSelection/rawCounts_modifiedNames/rawReads_combined.txt", header = T, sep = "\t") #not normalized
  
  #Convert gRNA into row names
  rownames(t01_allSamples_inputData) = t01_allSamples_inputData$gRNA
  t01_allSamples_inputData$gRNA = NULL
  
  #Select columns for analysis
  t01_allSamples_inputData = t01_allSamples_inputData %>%
    dplyr::select(one_of(c(SelectedCellLine_ungated, SelectedCellLine_testSample)))

  #Commentary
  print(paste("Library contains", nrow(t01_allSamples_inputData), "sgRNAs", "(including controls)", sep = " "))
  flush.console()
  
  #### Determine sequencing depth filter ###
    #Determine a sequencing depth filter. This helps me remove sgRNAs with low representation from the screen
    seqDepthFilter = t01_allSamples_inputData_rawCounts %>%
      gather("sample", "counts", 2:ncol(t01_allSamples_inputData_rawCounts))
    seqDepthFilter[is.na(seqDepthFilter)] = 0
    seqDepthFilter = seqDepthFilter %>%
      group_by(sample) %>%
      mutate(totalCounts = sum(counts)) %>%
      mutate(rpmFilter = 10 * (10^6 / totalCounts)) %>%
      ungroup(sample)
    
    seqDepthFilter = seqDepthFilter %>%
      dplyr::select(sample, rpmFilter) %>%
      distinct()
    
    seqDepthFilter_baseline = seqDepthFilter %>%
      filter(sample == SelectedCellLine_ungated)
    seqDepthFilter_baseline = seqDepthFilter_baseline$rpmFilter
    
    seqDepthFilter_testSample = seqDepthFilter %>%
      filter(sample == SelectedCellLine_testSample)
    seqDepthFilter_testSample = seqDepthFilter_testSample$rpmFilter

  ##################################
    
  #convert row names to column
  t01_allSamples_inputData$gRNA = rownames(t01_allSamples_inputData)
  rownames(t01_allSamples_inputData) = NULL
  
  #Add metadata to main file
  t01_allSamples_inputData = left_join(t01_allSamples_inputData, t01_metadata, by = "gRNA")
  
  #Detemrine the fold change of a sample compared to the baseline. Note that I add 0.01 reads to avoid 0 or inf values.
  t01_allSamples_inputData = t01_allSamples_inputData %>%
    mutate(FC = ((eval(as.symbol(SelectedCellLine_testSample)) + 0.01)/(eval(as.symbol(SelectedCellLine_ungated)) + 0.01))) #baseline is the melanoma population the day of the experiment without any selection
  
  #Calculate total number of guides per domain
  t01_allSamples_inputData = t01_allSamples_inputData %>%
    group_by(target) %>%
    mutate(numGuidesDesigned = n_distinct(gRNA)) %>%
    ungroup(target)

  #### Detemrine how many sgRNAs pass the sequencing depth filter ###
    tmp = t01_allSamples_inputData %>%
      filter(((eval(as.symbol(SelectedCellLine_testSample)) >= seqDepthFilter_testSample) | (eval(as.symbol(SelectedCellLine_ungated)) >= seqDepthFilter_baseline)))
    
    tmp = tmp %>%
      group_by(target) %>%
      mutate(numGuidesOverRPMfilter = n_distinct(gRNA)) %>%
      ungroup(target) %>%
      dplyr::select(target, numGuidesOverRPMfilter) %>%
      distinct()
    
    #Add information to each target about how many sgRNAs pass the sequencing depth filter.
    t01_allSamples_inputData = left_join(t01_allSamples_inputData, tmp, by = c("target"))

    
  ####################################################################
  
  ###NORMALIZATION###
  #We normalize data to the median effect observed by all sgRNAs in the screen. 
  
    #Calculate a value to use for normalization: median fold change across all non-targeting sgRNAs
    tmp = t01_allSamples_inputData %>%
      filter(target == "neg") %>%
      mutate(normalizationFactor = median(FC)) %>%
      dplyr::select(normalizationFactor) %>%
      distinct()
    
    #Commentary
    print(paste("Scaling factor = ", tmp$normalizationFactor, sep = " "))
    
    #Add the normalization value to the original dataset. 
    t01_allSamples_inputData = t01_allSamples_inputData %>%
      mutate(scalingFactor = tmp$normalizationFactor)
    
    #Normalize the fold change values and convert them to log2 values
    t01_allSamples_inputData = t01_allSamples_inputData %>%
      mutate(FC_normalized = FC / scalingFactor) %>%
      mutate(lFC = log2(FC_normalized))
  
    ####################################################################
    
    #Calculate normalized median FC after eliminating guides with low sequencing depth. We use this to obtain overall effect which we use for labeling in our plots.
    tmp = t01_allSamples_inputData %>%
      filter(((eval(as.symbol(SelectedCellLine_testSample)) >= seqDepthFilter_testSample) | (eval(as.symbol(SelectedCellLine_ungated)) >= seqDepthFilter_baseline))) %>%
      group_by(target) %>%
      mutate(medianlFC = median(lFC)) %>%
      mutate(sdlFC = sd(lFC)) %>%
      ungroup(target) %>%
      dplyr::select(target, medianlFC, sdlFC) %>%
      distinct()
    
    t01_allSamples_inputData = left_join(t01_allSamples_inputData, tmp, by = "target")
    
    #Calculate the number of guides over fold change filter
    tmp = t01_allSamples_inputData %>%
      filter(((eval(as.symbol(SelectedCellLine_testSample)) >= seqDepthFilter_testSample) | (eval(as.symbol(SelectedCellLine_ungated)) >= seqDepthFilter_baseline))) %>%
      filter(lFC >= 1 | lFC <= -1) %>%
      group_by(target) %>%
      mutate(nGuidesOverFC = n_distinct(gRNA)) %>%
      ungroup(target)
    tmp = tmp %>%
      dplyr::select(target, nGuidesOverFC) %>%
      distinct()
    
    t01_allSamples_inputData = left_join(t01_allSamples_inputData, tmp, by = "target")
    t01_allSamples_inputData$nGuidesOverFC[is.na(t01_allSamples_inputData$nGuidesOverFC)] = 0
    
    #Determine the ratio of guides that pass the fold change filter.
    t01_allSamples_inputData = t01_allSamples_inputData %>%
      mutate(fracGuidesOverFC = nGuidesOverFC/numGuidesOverRPMfilter)
 
   ###################################################################
  
  ### Directionality filter ###
    #Add directionality filter that eliminates targets whose guides are "hits" in both the enrichment side as well as in the depletion side
    tmp = t01_allSamples_inputData %>%
      filter(lFC >= 1 | lFC <= -1 ) %>%
      group_by(target) %>%
      mutate(PassDirectionalityFilter = ifelse((abs(max(lFC)) + abs(min(lFC)) == abs(max(lFC) + min(lFC))), T, F)) %>%
      ungroup(target) %>%
      dplyr::select(target, PassDirectionalityFilter) %>%
      distinct()
    # tmp$PassDirectionalityFilter[is.na(tmp$PassDirectionalityFilter)] = FALSE
    t01_allSamples_inputData <- left_join(t01_allSamples_inputData, tmp, by = "target")
  
    
    hits_066 = t01_allSamples_inputData %>%
      filter(fracGuidesOverFC >= 0.66) %>%
      filter(PassDirectionalityFilter == TRUE) %>%
      filter(!target %in% positiveControls)
    hits_066 = hits_066 %>%
      mutate(hitStatus_066 = "yes") %>%
      dplyr::select(target, medianlFC, sdlFC, hitStatus_066) %>%
      distinct()
  
  ###################################################################
  
  ### select list of hits to plot ###
  hits = hits_066 # 66% of sgRNAs  present in the screen must pass the effect size filter and the directionality filter
  
  #Save analyzed data
  fileName_allInfo = paste("03_extractedData/190122_targetedCrisprScreen_hitSelection/", currentSample, "analysis_", SelectedCellLine, "_allInfo.txt", sep = "")
  fileName_hitList = paste("03_extractedData/190122_targetedCrisprScreen_hitSelection/", currentSample, "analysis_", SelectedCellLine, "_hitList.txt", sep = "")
  write.table(t01_allSamples_inputData, file = fileName_allInfo, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(hits, file = fileName_hitList, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  

}
