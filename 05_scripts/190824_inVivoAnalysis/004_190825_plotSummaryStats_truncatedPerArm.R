#Make plots comparing KO to control until a specific timepoint

# # #For testing
# currentKO = "BRD2"

#Select minimum number of mice for plotting both controls and KOs (per treatment arm)
min_nMice = 3
max_time = 90

#Load data
allData = read.table("03_extractedData/190824_inVivoAnalysis/changeInTumorVolume.txt", header = T)

#Add mean day of measurement to the treatment groups
tmp = allData %>%
  dplyr::select(measurementGroup, treatmentDuration) %>%
  distinct() %>%
  group_by(measurementGroup) %>%
  mutate(treatmentDurationMean = mean(treatmentDuration)) %>%
  ungroup(measurementGroup)

#Add mean time of measurement group
allData = left_join(allData, tmp, by = c("measurementGroup", "treatmentDuration"))

#Calculate summary statistics
allData = allData %>%
  group_by(KO, treatmentArm, measurementGroup) %>%
  mutate(nMice = n_distinct(mouseID)) %>%
  mutate(meanVolume = mean(tumorVolume)) %>%
  mutate(sdVolume = sd(tumorVolume)) %>%
  mutate(seVolume = sdVolume / sqrt(nMice)) %>%
  ungroup(KO, treatmentArm, measurementGroup)

#Determine the different KOs
listOfKO = c("BRD2", "LATS2", "DOT1L")
listOfTreatmentArms = c("untreated", "treated")

for(currentKO in listOfKO) {
  
  ##Detemrine last timepoint to plot
  
  #Select time point to stop each comparison
  tmp_allData = allData %>%
    filter(KO == "Neg" | KO == currentKO)
  
  tmp_allData_endpoint = tmp_allData %>%
    dplyr::select(KO, treatmentArm, treatmentDurationMean, nMice) %>%
    distinct()
  
  tmp_allData_endpoint = tmp_allData_endpoint %>%
    filter(nMice >= min_nMice) %>%
    filter(treatmentDurationMean <= max_time) %>% #this line stops the plot at d30 b/c BRD2 mice that are responding very well are still not at d30.
    group_by(KO, treatmentArm) %>%
    mutate(endPoint = max(treatmentDurationMean)) %>%
    ungroup(KO, treatmentArm)
  
  tmp_allData_endpoint = tmp_allData_endpoint %>%
    dplyr::select(KO, treatmentArm, endPoint) %>%
    distinct()
  
  tmp_allData_endpoint = tmp_allData_endpoint %>%
    group_by(treatmentArm) %>%
    mutate(endPoint_joint = min(endPoint)) %>%
    ungroup(treatmentArm) %>%
    dplyr::select(-endPoint, -KO) %>%
    distinct()
  
  #Trucante data at the point and plot
  allData_truncated = left_join(tmp_allData, tmp_allData_endpoint, by = c("treatmentArm"))
  allData_truncated = allData_truncated %>%
    filter(treatmentDurationMean <= endPoint_joint + 1) #I add one because sometimes the groups are a day apart in measurement
  
  growthPlot_onDrug = ggplot(data = allData_truncated, aes(x = treatmentDurationMean, y = meanVolume, color = KO)) +
    geom_point(size = 3) +
    geom_line() +
    geom_errorbar(aes(ymin = (meanVolume - seVolume), ymax = (meanVolume + seVolume), width = 0.5)) +
    facet_grid(treatmentArm~.) +
    xlab("days on treatment") +
    ylab("tumor volume") +
    xlim(0,40) +
    ylim(0,2500) +
    ggtitle(paste(currentKO, "_tumorVolumes", sep = "")) +
    #geom_hline(yintercept = 0) +
    theme_classic()
  growthPlot_onDrug
  plotName = paste("04_plots/190824_inVivoAnalysis/",currentKO, "_lineplot_tuncated.PDF", sep = "")
  ggsave(plotName)
  
}

