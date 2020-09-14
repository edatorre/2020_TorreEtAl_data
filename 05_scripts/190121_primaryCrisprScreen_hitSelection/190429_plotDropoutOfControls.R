# #The goal of this script is to plot the effect of negative and positive controls over time in my screens

#Read in data
epigeneticData = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_epigenetic_allInfo.txt", header = T, sep = "\t")

#Make list of positive controls
positiveControls = c("PCNA", "POLR2A", "POLR2D", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")

#Obtain normalizationFactor
normalizationFactor = epigeneticData %>%
  filter(gene == "neg") 

normalizationFactor = normalizationFactor %>%
  summarize(normFactor = mean(DropOut_FC))
  
#Isolate controls
epigeneticData_negControls = epigeneticData %>%
  filter(gene == "neg") %>%
  mutate(sampleType = "negativeControl")

epigeneticData_PosControls = epigeneticData %>%
  filter(gene %in% positiveControls) %>%
  mutate(sampleType = "positiveControl")

medianFC_negativeControls = epigeneticData_negControls %>%
  summarize(medianOfNegControls = median(DropOut_FC))

#Merge data
epigeneticData_controls = bind_rows(epigeneticData_negControls, epigeneticData_PosControls)

#Remove NA rows
epigeneticData_controls = epigeneticData_controls[is.finite(epigeneticData_controls$DropOut_FC), ]

#Normalize FC values to the negative controls
epigeneticData_controls = epigeneticData_controls %>%
  dplyr::select(gRNA, DropOut_FC, sampleType)
epigeneticData_controls$scalingFactor = medianFC_negativeControls$medianOfNegControls
epigeneticData_controls = epigeneticData_controls %>%
  mutate(DropOut_FC_norm = DropOut_FC / scalingFactor)
epigeneticData_controls = epigeneticData_controls %>%
  mutate(DropOut_lFC_norm = log2(DropOut_FC_norm + 0.01))

epigeneticData_controls = epigeneticData_controls %>%
  group_by(sampleType) %>%
  mutate(medianlFC_norm = median(DropOut_lFC_norm)) %>%
  ungroup(sampleType) %>%
  arrange(sampleType, gRNA) %>%
  mutate(plottingID = 1:nrow(epigeneticData_controls))

medianFC = epigeneticData_controls %>%
  dplyr::select(sampleType, medianlFC_norm) %>%
  distinct()

#Create plots
controlsPlot = ggplot(data = epigeneticData_controls, aes(x = plottingID, y = DropOut_lFC_norm, fill = sampleType)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = medianFC$medianlFC_norm, linetype = "dashed") +
  theme_classic() +
  ylab("Log2 Fold change over time without selection")
controlsPlot
ggsave("04_plots/190121_primaryCrisprScreen_hitSelection/epigeneticLibrary_controls_dropout.PDF") 


##########

#Read in data
kinaseData = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_kinase_allInfo.txt", header = T, sep = "\t")

#Make list of positive controls
positiveControls = c("PCNA", "POLR2A", "POLR2D", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")

#Obtain normalizationFactor
normalizationFactor = kinaseData %>%
  filter(gene == "neg") 

normalizationFactor = normalizationFactor %>%
  summarize(normFactor = mean(DropOut_FC))

#Isolate controls
kinaseData_negControls = kinaseData %>%
  filter(gene == "neg") %>%
  mutate(sampleType = "negativeControl")

kinaseData_PosControls = kinaseData %>%
  filter(gene %in% positiveControls) %>%
  mutate(sampleType = "positiveControl")

medianFC_negativeControls = kinaseData_negControls %>%
  summarize(medianOfNegControls = median(DropOut_FC))

#Merge data
kinaseData_controls = bind_rows(kinaseData_negControls, kinaseData_PosControls)

#Remove NA rows
kinaseData_controls = kinaseData_controls[is.finite(kinaseData_controls$DropOut_FC), ]

#Normalize FC values to the negative controls
kinaseData_controls = kinaseData_controls %>%
  dplyr::select(gRNA, DropOut_FC, sampleType)
kinaseData_controls$scalingFactor = medianFC_negativeControls$medianOfNegControls
kinaseData_controls = kinaseData_controls %>%
  mutate(DropOut_FC_norm = DropOut_FC / scalingFactor)
kinaseData_controls = kinaseData_controls %>%
  mutate(DropOut_lFC_norm = log2(DropOut_FC_norm + 0.01))

kinaseData_controls = kinaseData_controls %>%
  group_by(sampleType) %>%
  mutate(medianlFC_norm = median(DropOut_lFC_norm)) %>%
  ungroup(sampleType) %>%
  arrange(sampleType, gRNA) %>%
  mutate(plottingID = 1:nrow(kinaseData_controls))

medianFC = kinaseData_controls %>%
  dplyr::select(sampleType, medianlFC_norm) %>%
  distinct()

#Create plots
controlsPlot = ggplot(data = kinaseData_controls, aes(x = plottingID, y = DropOut_lFC_norm, fill = sampleType)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = medianFC$medianlFC_norm, linetype = "dashed") +
  theme_classic() +
  ylab("Log2 Fold change over time without selection")
controlsPlot
ggsave("04_plots/190121_primaryCrisprScreen_hitSelection/kinaseLibrary_controls_dropout.PDF") 


#######

#Read in data
tfData = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/DPanalysis_tf_allInfo.txt", header = T, sep = "\t")

#Make list of positive controls
positiveControls = c("PCNA", "POLR2A", "POLR2D", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")

#Obtain normalizationFactor
normalizationFactor = tfData %>%
  filter(gene == "neg") 

normalizationFactor = normalizationFactor %>%
  summarize(normFactor = mean(DropOut_FC))

#Isolate controls
tfData_negControls = tfData %>%
  filter(gene == "neg") %>%
  mutate(sampleType = "negativeControl")

tfData_PosControls = tfData %>%
  filter(gene %in% positiveControls) %>%
  mutate(sampleType = "positiveControl")

medianFC_negativeControls = tfData_negControls %>%
  summarize(medianOfNegControls = median(DropOut_FC))

#Merge data
tfData_controls = bind_rows(tfData_negControls, tfData_PosControls)

#Remove NA rows
tfData_controls = tfData_controls[is.finite(tfData_controls$DropOut_FC), ]

#Normalize FC values to the negative controls
tfData_controls = tfData_controls %>%
  dplyr::select(gRNA, DropOut_FC, sampleType)
tfData_controls$scalingFactor = medianFC_negativeControls$medianOfNegControls
tfData_controls = tfData_controls %>%
  mutate(DropOut_FC_norm = DropOut_FC / scalingFactor)
tfData_controls = tfData_controls %>%
  mutate(DropOut_lFC_norm = log2(DropOut_FC_norm + 0.01))

tfData_controls = tfData_controls %>%
  group_by(sampleType) %>%
  mutate(medianlFC_norm = median(DropOut_lFC_norm)) %>%
  ungroup(sampleType) %>%
  arrange(sampleType, gRNA) %>%
  mutate(plottingID = 1:nrow(tfData_controls))

medianlFC = tfData_controls %>%
  dplyr::select(sampleType, medianlFC_norm) %>%
  distinct()

#Create plots
controlsPlot = ggplot(data = tfData_controls, aes(x = plottingID, y = DropOut_lFC_norm, fill = sampleType)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = medianlFC$medianlFC_norm, linetype = "dashed") +
  theme_classic() +
  ylab("Log2 Fold change over time without selection")
controlsPlot
ggsave("04_plots/190121_primaryCrisprScreen_hitSelection/tfLibrary_controls_dropout.PDF") 
