#This function selects hits from the screen and plots the results. 

# library(tidyverse)
# library(RColorBrewer)

selectHits_primaryScreen = function(CRISPR_library, BiologicalSampleName) {

  # #For testing
  # CRISPR_library = "TF"
  # BiologicalSampleName = "DP"
  
  #Determine variables
  currentLibrary = CRISPR_library
  currentSample = BiologicalSampleName
  currentSample_FC = paste(currentSample, "_FC", sep = "")
  currentSample_lFC = paste(currentSample, "_lFC", sep = "")
  
  #Define the positive controls
  positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
  negativeControls = c("neg", "WT")
  
  #Read in metadata that contains information which GO term a given KO is part of. I will use this for plotting later on. 
  metadataFileName = paste("03_extractedData/190121_primaryCrisprScreen_hitSelection/", currentLibrary, "Targets_assignedToGOs.txt", sep = "")
  t01_metadata = read.table(metadataFileName, sep = "\t", header = T)
  
  #Read screen data. These are the rpm-normalized read counts for each guide and well as the unnormalized values.
  t01_allSamples_inputData = read.table("01_rawData/190121_primaryCrisprScreen_hitSelection/rawReads/CRISPRscreens_combined_rpmNormalized.txt", header = T, sep = "\t") #normalized
  t01_allSamples_inputData_rawCounts = read.table("01_rawData/190121_primaryCrisprScreen_hitSelection/rawReads/CRISPRscreens_combined.txt", header = T, sep = "\t") #not normalized
  
  #Select library for analysis
  t01_allSamples_inputData = t01_allSamples_inputData %>%
    filter(library == currentLibrary)
  t01_allSamples_inputData_rawCounts = t01_allSamples_inputData_rawCounts %>%
    filter(library == currentLibrary)
  
  #Commentary
  numGuidesInLib = nrow(t01_allSamples_inputData)
  print(paste("Library contains", numGuidesInLib, "sgRNAs", "(including controls)", sep = " "))
  flush.console()
  
  #Obtain total number of domains targeted
  tmp = t01_allSamples_inputData %>%
    dplyr::select(domain) %>%
    distinct()
  
  #Commentary
  print(paste("It targets", nrow(tmp), "different proteins", sep = " "))
  flush.console()
  
  #### Determine sequencing depth filter ###
    #Determine a sequencing depth filter. This helps me remove sgRNAs with low representation from the screen
    seqDepthFilter = t01_allSamples_inputData_rawCounts %>%
      gather("sample", "counts", 6:11)
    seqDepthFilter[is.na(seqDepthFilter)] = 0
    seqDepthFilter = seqDepthFilter %>%
      filter(sample == currentSample | sample == "ungated") ####Here I am selecting the filter based on the ungated sample and the test sample
    seqDepthFilter = seqDepthFilter %>%
      group_by(sample) %>%
      mutate(totalCounts = sum(counts)) %>%
      ungroup(sample)
    
    totalReads_control = seqDepthFilter %>%
      filter(sample == "ungated") %>%
      dplyr::select(sample, totalCounts) %>%
      distinct()
    totalReads_sample = seqDepthFilter %>%
      filter(sample == `currentSample`) %>%
      dplyr::select(sample, totalCounts) %>%
      distinct()
    
    seqDepthFilter_rpm_control = 10 * (10^6 / totalReads_control$totalCounts)
    seqDepthFilter_rpm_sample = 10 * (10^6 / totalReads_sample$totalCounts)
    
    #Add sequencing depth filters to main data
    t01_allSamples_inputData = t01_allSamples_inputData %>%
      mutate(seqDepthFilter_unagated = seqDepthFilter_rpm_control) %>%
      mutate(seqDepthFilter_sample = seqDepthFilter_rpm_sample)
  ##################################
  
  #Add metadata to main file
  t01_allSamples_inputData = left_join(t01_allSamples_inputData, t01_metadata, by = "gene")
  
  #Calculate total number of guides per domain
  t01_allSamples_inputData = t01_allSamples_inputData %>%
    group_by(domain) %>%
    mutate(numGuidesDesigned = n_distinct(gRNA)) %>%
    ungroup(domain)
  
  t01_allSamples_inputData = t01_allSamples_inputData %>%
    mutate(numGuidesDesgined_Filter = ifelse(numGuidesDesigned >= 4, "yes", "no"))
  
  #comment on number of domains with less than 4 sgRNAs
  tmp = t01_allSamples_inputData %>%
    filter(numGuidesDesgined_Filter == "no") %>%
    dplyr::select(domain) %>%
    distinct()
  
  print(paste(nrow(tmp), "domains have less than 4 sgRNAs designed. Those are removed from the analysis"))
    
  #Identify sgRNAs with infinite fold changes I will use this for plotting later. Add this info to main file
  infiniteFC_ungated = t01_allSamples_inputData %>%
    filter(ungated == 0 & eval(as.symbol(currentSample)) != 0) %>%
    mutate(infiniteFC = "yes")
  infiniteFC_ungated = infiniteFC_ungated %>%
    dplyr::select(gRNA, infiniteFC)
  
  infiniteFC_sample = t01_allSamples_inputData %>%
    filter(eval(as.symbol(currentSample)) == 0 & ungated != 0) %>%
    mutate(infiniteFC = "yes")
  infiniteFC_sample = infiniteFC_sample %>%
    dplyr::select(gRNA, infiniteFC)
  
  infiniteFC_sgRNA = bind_rows(infiniteFC_ungated, infiniteFC_sample)
  
  t01_allSamples_inputData = left_join(t01_allSamples_inputData, infiniteFC_sgRNA, by = "gRNA")
  t01_allSamples_inputData$infiniteFC[is.na(t01_allSamples_inputData$infiniteFC)] = "no"
  
  #Commentary
  print(paste(nrow(infiniteFC_sgRNA), "sgRNAs have an infinite fold change in the screen", sep = " "))
  flush.console()
  
  #### Detemrine how many sgRNAs pass the sequencing depth filter ####
  tmp = t01_allSamples_inputData_rawCounts %>%
    filter(ungated >= t01_allSamples_inputData$seqDepthFilter_unagated | (eval(as.symbol(currentSample))) >= t01_allSamples_inputData$seqDepthFilter_sample) %>%
    group_by(domain) %>%
    mutate(numGuidesOverRPMfilter = n_distinct(gRNA)) %>%
    ungroup(domain) %>%
    dplyr::select(domain, numGuidesOverRPMfilter) %>%
    distinct()
  
  tmp2 = t01_allSamples_inputData_rawCounts %>%
    filter(ungated < t01_allSamples_inputData$seqDepthFilter_unagated & (eval(as.symbol(currentSample))) < t01_allSamples_inputData$seqDepthFilter_sample) %>%
    mutate(representation = "low") %>%
    dplyr::select(gRNA, representation)
  
  #Commentary
  print(paste(nrow(tmp2), " sgRNAs fail the sequencing depth filter"))
  flush.console()
  
  #add adata to main file
  t01_allSamples_inputData = left_join(t01_allSamples_inputData, tmp, by = "domain")
  t01_allSamples_inputData = left_join(t01_allSamples_inputData, tmp2, by = "gRNA")
  t01_allSamples_inputData$representation[is.na(t01_allSamples_inputData$representation)] = "normal"

  ####################################################################  
  ###FOLD CHANGE###
    
  #Detemrine the fold change of a sample compared to the baseline. Note that I add 0.01 reads to avoid 0 or inf values.
  t01_allSamples_inputData = t01_allSamples_inputData %>%
    mutate(DropOut_FC = ((ungated + 0.01)/(baseline + 0.01))) %>% #this tells me if a sgRNA lost/gained representation without NGFR/EGFR of vemurafenib resistance
    mutate(FC = ((eval(as.symbol(currentSample)) + 0.01)/(ungated + 0.01))) #baseline is the melanoma population the day of the experiment without any selection

  t01_allSamples_inputData$FC[is.na(t01_allSamples_inputData$FC)] = 0 #without this I cannot obtain a median FC for normalization. 
  ####################################################################
  
  ###NORMALIZATION###
  #We normalize data to the median effect observed by all sgRNAs in the screen. 
  
    #Calculate a value to use for normalization: median fold change across all sgRNAs
    tmp = t01_allSamples_inputData %>%
      group_by(library) %>%
      mutate(normalizationFactor = median(FC)) %>%
      ungroup(library) %>%
      dplyr::select(library, normalizationFactor) %>%
      distinct()
    
    #Commentary
    print(paste("Scaling factor = ", tmp$normalizationFactor, sep = " "))
    
    #Add the normalization value to the original dataset. 
    t01_allSamples_inputData = left_join(t01_allSamples_inputData, tmp, by = "library")
    
    #Normalize the fold change values and convert them to log2 values
    t01_allSamples_inputData = t01_allSamples_inputData %>%
      mutate(FC_normalized = FC / normalizationFactor) %>%
      mutate(lFC = log2(FC_normalized))
  ####################################################################

    #Calculate normalized median FC after eliminating guides with low sequencing depth. We use this to obtain overall effect which we use for labeling in our plots.
    tmp = t01_allSamples_inputData %>%
      filter(representation != "low") %>%
      group_by(domain) %>%
      mutate(medianlFC = median(lFC)) %>%
      ungroup(domain) %>%
      dplyr::select(domain, library, medianlFC) %>%
      distinct()
    
    t01_allSamples_inputData = left_join(t01_allSamples_inputData, tmp, by = c("domain", "library"))
  ####################################################################
  
  ### Effect-size filters ###
  
    #Determine number of guides above or below a two-fold change
    tmp <- t01_allSamples_inputData %>%
      filter(lFC >= 1 | lFC <= -1 ) %>%
      filter(representation != "low") %>%
      group_by(domain, library) %>%
      mutate(nGuidesOverFC = n_distinct(gRNA)) %>%
      ungroup(domain, library) %>%
      dplyr::select(domain, nGuidesOverFC) %>%
      distinct()
    
    #Add the number of guides passing the fold change filter to the main sample file
    t01_allSamples_inputData <- left_join(t01_allSamples_inputData, tmp, by = "domain")
    t01_allSamples_inputData$nGuidesOverFC[is.na(t01_allSamples_inputData$nGuidesOverFC)] = 0 #convert NAs to 0.
    
    #Determine the ratio of guides that pass the fold change filter.
    t01_allSamples_inputData = t01_allSamples_inputData %>%
      mutate(fracGuidesOverFC = nGuidesOverFC/numGuidesOverRPMfilter)
  ###################################################################
  
  ### Directionality filter ###
    #Add directionality filter that eliminates targets whose guides are "hits" in both the enrichment side as well as in the depletion side
    tmp = t01_allSamples_inputData %>%
      filter(lFC >= 1 | lFC <= -1 ) %>%
      filter(representation != "low") %>%
      group_by(domain) %>%
      mutate(PassDirectionalityFilter = ifelse((abs(max(lFC)) + abs(min(lFC)) == abs(max(lFC) + min(lFC))), T, F)) %>%
      ungroup(domain) %>%
      dplyr::select(domain, PassDirectionalityFilter) %>%
      distinct()
    tmp$PassDirectionalityFilter[is.na(tmp$PassDirectionalityFilter)] = FALSE
    t01_allSamples_inputData <- left_join(t01_allSamples_inputData, tmp, by = "domain")
  
  ### Select hits ###
    #Keep only genes with > x guides passing filter in the same direction
    hits_050 = t01_allSamples_inputData %>%
      filter(fracGuidesOverFC >= 0.50) %>%
      filter(PassDirectionalityFilter == TRUE) %>%
      filter(numGuidesDesgined_Filter == "yes") %>%
      filter(!gene %in% negativeControls) %>%
      filter(!gene %in% positiveControls)
    hits_050 = hits_050 %>%
      mutate(hitStatus_050 = "yes") %>%
      dplyr::select(domain, library, medianlFC, hitStatus_050) %>%
      distinct()
    
    hits_066 = t01_allSamples_inputData %>%
      filter(fracGuidesOverFC >= 0.66) %>%
      filter(PassDirectionalityFilter == TRUE) %>%
      filter(numGuidesDesgined_Filter == "yes") %>%
      filter(!gene %in% negativeControls) %>%
      filter(!gene %in% positiveControls)
    hits_066 = hits_066 %>%
      mutate(hitStatus_066 = "yes") %>%
      dplyr::select(domain, library, medianlFC, hitStatus_066) %>%
      distinct()
    
    hits_075 = t01_allSamples_inputData %>%
      filter(fracGuidesOverFC >= 0.75) %>%
      filter(PassDirectionalityFilter == TRUE) %>%
      filter(numGuidesDesgined_Filter == "yes") %>%
      filter(!gene %in% negativeControls) %>%
      filter(!gene %in% positiveControls)
    hits_075 = hits_075 %>%
      mutate(hitStatus_075 = "yes") %>%
      dplyr::select(domain, library, medianlFC, hitStatus_075) %>%
      distinct()
    
    hits_090 = t01_allSamples_inputData %>%
      filter(fracGuidesOverFC >= 0.90) %>%
      filter(PassDirectionalityFilter == TRUE) %>%
      filter(numGuidesDesgined_Filter == "yes") %>%
      filter(!gene %in% negativeControls) %>%
      filter(!gene %in% positiveControls)
    hits_090 = hits_090 %>%
      mutate(hitStatus_090 = "yes") %>%
      dplyr::select(domain, library, medianlFC, hitStatus_090) %>%
      distinct()
  ###################################################################
  
  ### select list of hits to plot ###
  hits = hits_066 # 66% of sgRNAs  present in the screen must pass the effect size filter and the directionality filter
  hits_up = hits %>%
    filter(medianlFC > 0)
  hits_down = hits %>%
    filter(medianlFC < 0)
    
    #Commentary
    print(paste("There are ", nrow(hits), " hits (by domain) in the ", CRISPR_library," library" , sep = ""))
    print(paste(nrow(hits_up), " increase the ", BiologicalSampleName, " phenotype", sep = ""))
    print(paste(nrow(hits_down), " decrease the ", BiologicalSampleName, " phenotype", sep = ""))
  
  #Save analyzed data
  fileName_allInfo = paste("03_extractedData/190121_primaryCrisprScreen_hitSelection/", currentSample, "analysis_", currentLibrary, "_allInfo.txt", sep = "")
  fileName_hitList = paste("03_extractedData/190121_primaryCrisprScreen_hitSelection/", currentSample, "analysis_", currentLibrary, "_hitList.txt", sep = "")
  write.table(t01_allSamples_inputData, file = fileName_allInfo, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(hits, file = fileName_hitList, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  ####PLOT RESULTS####
  
    #Make a table that contians only the data that we will plot
    #Work on negative and positive controls
    plottingTable_1_negativeControls = t01_allSamples_inputData %>%
      filter(domain == "neg") %>%
      arrange(ID, library, domain)
    
    plottingTable_1_negativeControls$Description = gsub("other", "negativeControl", plottingTable_1_negativeControls$Description)
    
    print(paste("Plot contains", nrow(plottingTable_1_negativeControls), "non-targeting sgRNAs as negative controls", sep = " "))
    flush.console
    
    plottingTable_1_positiveControls = t01_allSamples_inputData %>%
      filter(gene %in% positiveControls) %>%
      arrange(ID, library, domain)
    
    plottingTable_1_positiveControls$Description = gsub("other", "positiveControl", plottingTable_1_positiveControls$Description)
    
    print(paste("Plot contains", nrow(plottingTable_1_positiveControls), "sgRNAs as cell viability controls", sep = " "))
    flush.console
    
    #Work on the rest of the data
    plottingTable_1 = t01_allSamples_inputData %>%
      filter(library == currentLibrary) %>%
      filter(domain != "neg") %>% #remove positive and negative controls from the dataset
      filter(!gene %in% positiveControls) %>%
      arrange(ID, library, domain) 
    
    plottingTable_1 = plottingTable_1 %>%
      distinct()
    
    plottingTable_1 = plottingTable_1 %>%
      mutate(plottingID = 1:nrow(plottingTable_1)) #add a number that will allow to plots the sgRNAs in a specific order
    
    firstPlotID_positiveControls = max(plottingTable_1$plottingID) + 300 #these number simply add spacing between the test samples and the controls.
    lastPlotID_positiveControls = firstPlotID_positiveControls + nrow(plottingTable_1_positiveControls)
    firstPlotID_negativeControls = lastPlotID_positiveControls + 300
    lastPlotID_negativeControls = firstPlotID_negativeControls + nrow(plottingTable_1_negativeControls)
    
    plottingTable_1_positiveControls = plottingTable_1_positiveControls %>%
      mutate(plottingID = firstPlotID_positiveControls:(lastPlotID_positiveControls - 1))
    
    plottingTable_1_negativeControls = plottingTable_1_negativeControls %>%
      mutate(plottingID = firstPlotID_negativeControls:(lastPlotID_negativeControls - 1))
    
    plottingTable_1 = bind_rows(plottingTable_1, plottingTable_1_positiveControls, plottingTable_1_negativeControls)
  
    #Add label describing whether a target is a hit or not
    plottingTable_1 = plottingTable_1 %>%
      mutate(geneHighlight = ifelse(domain %in% hits$domain, TRUE, FALSE))
    
    #Make table with only the names of the hits
    plottingTable_1_namesOfHits = plottingTable_1 %>%
      filter(domain %in% hits$domain)
    
    #Keep the maximum value of the hits. That vaue will determine the coordinate of the text.
    plottingTable_1_namesOfHits = plottingTable_1_namesOfHits %>%  
      group_by(domain) %>%
      mutate(maxAbsVal = max(abs(lFC))) %>%
      ungroup(domain) %>%
      filter(maxAbsVal == abs(lFC)) %>%
      dplyr::select(domain, lFC, Description, plottingID, geneHighlight) %>%
      distinct()
    
    #save this table to keep the order of the hits and the groups
    write.table(plottingTable_1_namesOfHits, file = paste("03_extractedData/190121_primaryCrisprScreen_hitSelection/", currentSample, "_", currentLibrary, "Library_groupOrders.txt", sep = ""), 
                col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

    #Keep the maximum value of the controls. That vaue will determine the coordinate of the text.
    plottingTable_1_postivieControls_names = plottingTable_1_positiveControls %>%  
      group_by(gene) %>%
      mutate(maxAbsVal = max(abs(lFC))) %>%
      ungroup(gene) %>%
      filter(maxAbsVal == abs(lFC)) %>%
      dplyr::select(domain, lFC, Description, plottingID) %>%
      distinct()
    
    #Determine the coordinate where a given category of targets end
    classLimits = plottingTable_1 %>%
      dplyr::select(Description, plottingID) %>%
      group_by(Description) %>%
      summarize(maxPLottingID = max(plottingID))
    
    classLimits = classLimits$maxPLottingID
    
    #Separate hits from non-hits to better color code the screen.
    plottingTable_1_hits = plottingTable_1 %>%
      filter(domain %in% hits$domain)
    
    plottingTable_1_nothits = plottingTable_1 %>%
      filter(!domain %in% hits$domain)

    
    #Plot the results: main plot
    plot_dotPlot <- ggplot() +
      geom_point(data = plottingTable_1_hits, shape = 20, aes(x = -plottingID, y = lFC, colour = Description, shape = infiniteFC), stroke = 0.3, size = 2) +
      geom_point(data = plottingTable_1_nothits, shape = 20, aes(x = -plottingID, y = lFC, shape = infiniteFC), colour = "gray", stroke = 0.2, size = 0.3) +
      geom_text(data = plottingTable_1_namesOfHits, aes(x = -plottingID, y = lFC, label = domain), size = 2) +
      #geom_text_repel(data = plottingTable_1_postivieControls_names, aes(x = -plottingID, y = lFC, label = domain), size = 2) +
      scale_shape_manual(values = c("yes" = 4, "no" = 16)) +
      theme_classic() +
      theme(legend.position="none") +
      scale_y_continuous(breaks = round(seq(-20, 10, by = 1),1)) +
      geom_hline(yintercept = c(1, -1), linetype="dashed") +
      geom_vline(xintercept = -1 * classLimits) +
      coord_flip()
    plot(plot_dotPlot)
    
    #Make plot name and save the plot
    plotName = paste("04_plots/190121_primaryCrisprScreen_hitSelection/", currentSample, "_", currentLibrary, "Library_pointPlot.PDF", sep = "")
    ggsave(plotName, plot = last_plot(), width = 3, height = (numGuidesInLib / 500))
    
    #Plot the results: highlights guides with infinite fold change
    plot_dotPlot_2 <- ggplot() +
      geom_point(data = plottingTable_1_hits, shape = 20, aes(x = -plottingID, y = lFC, colour = infiniteFC), stroke = 0.3, size = 2) +
      geom_point(data = plottingTable_1_nothits, shape = 20, aes(x = -plottingID, y = lFC, colour = infiniteFC), stroke = 0.2, size = 0.3) +
      geom_text(data = plottingTable_1_namesOfHits, aes(x = -plottingID, y = lFC, label = domain), size = 2) +
      #geom_text_repel(data = plottingTable_1_postivieControls_names, aes(x = -plottingID, y = lFC, label = domain), size = 2) +
      scale_colour_manual(values = c("yes" = "red", "no" = "grey")) +
      theme_classic() +
      theme(legend.position="none") +
      scale_y_continuous(breaks = round(seq(-20, 10, by = 1),1)) +
      geom_hline(yintercept = c(1, -1), linetype="dashed") +
      geom_vline(xintercept = -1 * classLimits) +
      coord_flip()
    plot(plot_dotPlot_2)
    
    #Make plot name and save the plot
    plotName = paste("04_plots/190121_primaryCrisprScreen_hitSelection/", currentSample, "_", currentLibrary, "Library_pointPlot_highlitghInfiniteFC.PDF", sep = "")
    ggsave(plotName, plot = last_plot(), width = 3, height = (numGuidesInLib / 500))
  
  ### TIERS ####
    ###organize hits by Tiers depending on how many sgRNAs pass the filters
    #Tiers: 75% sgRNAs pass filter -> Tier1, 66% -> Tier2, 50% -> Tier3, less than 50% -> Tier4
    
    #Limit dataset to useful fields
    slidingScaleTable = t01_allSamples_inputData %>%
      dplyr::select(gRNA, domain, library) %>%
      distinct()
    
    #Add data specifying if a target is a hit in a given tier
    slidingScaleTable = left_join(slidingScaleTable, hits_050, by = c("domain", "library"))
    slidingScaleTable = left_join(slidingScaleTable, hits_066, by = c("domain", "library", "medianlFC"))
    slidingScaleTable = left_join(slidingScaleTable, hits_075, by = c("domain", "library", "medianlFC"))
    
    slidingScaleTable$medianlFC = NULL
    slidingScaleTable[is.na(slidingScaleTable)] = "No" #Convert NA's to NO. 
    
    slidingScaleTable = slidingScaleTable %>% #Organize table to tiers. 
      arrange(desc(hitStatus_075), desc(hitStatus_066), desc(hitStatus_050))
    
    #save the table
    slidingScaleTableName = paste("03_extractedData/190121_primaryCrisprScreen_hitSelection/", currentSample, "_", currentLibrary,"Library_slidingScaleTable.txt", sep = "")
    write.table(slidingScaleTable, file = slidingScaleTableName, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    
}
