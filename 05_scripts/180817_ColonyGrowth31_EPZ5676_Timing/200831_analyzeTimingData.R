# #The goal of this script is to process colony growth data from EPX5676 timing experiments.

#Read in platemaps
plateMap = read.table("02_metadata/180817_ColonyGrowth31_EPZ5676_Timing/plateMap.txt", header = T)

#Read in summarize plate data
baseline = read.csv("03_extractedData/180817_ColonyGrowth31_EPZ5676_Timing/plateSummaries/baseline_plateSummary.csv")
plateA = read.csv("03_extractedData/180817_ColonyGrowth31_EPZ5676_Timing/plateSummaries/plateA_plateSummary.csv")
plateB = read.csv("03_extractedData/180817_ColonyGrowth31_EPZ5676_Timing/plateSummaries/plateB_plateSummary.csv")
plateC = read.csv("03_extractedData/180817_ColonyGrowth31_EPZ5676_Timing/plateSummaries/plateC_plateSummary.csv")

#Add plate info to each summary file
baseline = baseline %>%
  mutate(PlateName = "baseline")
plateA = plateA %>%
  mutate(PlateName = "plateA")
plateB = plateB %>%
  mutate(PlateName = "plateB")
plateC = plateC %>%
  mutate(PlateName = "plateC")

#Select relevant baseline data and rename columns
baseline = baseline %>%
  dplyr::select(well_number, num_cells)
colnames(baseline)[2] = "baseline_num_cells"

#Combine data files
allData = bind_rows(plateA, plateB, plateC)

#Add baseline data
allData = left_join(allData, baseline, by = "well_number")

#Add metadata
allData = inner_join(plateMap, allData, by = c("PlateName", "well_number"))

#Normalize values to the total number of cell plates in the well.
allData = allData %>%
  mutate(Rcells_norm = num_cells / baseline_num_cells) %>%
  mutate(RcellsOutsideColonies_norm = cells_outside_colonies / baseline_num_cells) %>%
  mutate(Rcolonies_norm = num_colonies * 1000 / baseline_num_cells) 

#Obtain values for control treatment
controlData_Rcolonies = allData %>%
  filter(SampleDescriptions == "DMSO_PLX4032") %>%
  summarize(mean(Rcolonies_norm)) 

#Obtain values for control treatment
controlData_RcellsOutsideColonies = allData %>%
  filter(SampleDescriptions == "DMSO_PLX4032") %>%
  summarize(mean(RcellsOutsideColonies_norm)) 

#Obtain fold change over control for the number of resistance colonies
allData = allData %>%
  mutate(Rcolonies_norm_fc = Rcolonies_norm / controlData_Rcolonies$`mean(Rcolonies_norm)`) %>%
  mutate(RcellsOutsideColonies_norm_fc = RcellsOutsideColonies_norm / controlData_RcellsOutsideColonies$`mean(RcellsOutsideColonies_norm)`)

#Obtain N per group
allData = allData %>%
  group_by(SampleDescriptions) %>%
  mutate(nSamples = n_distinct(PlateName)) %>%
  ungroup(SampleDescriptions)

#Generate means and SD for the different metrics
allData = allData %>%
  group_by(SampleDescriptions) %>%
  mutate(meanRcolonies = mean(Rcolonies_norm_fc)) %>%
  mutate(sdRcolonies = sd(Rcolonies_norm_fc)) %>%
  mutate(seRcolonies = sdRcolonies / sqrt(nSamples)) %>%
  mutate(meanRcellsOutsideColonies = mean(RcellsOutsideColonies_norm_fc)) %>%
  mutate(sdRcellsOutsideColonies = sd(RcellsOutsideColonies_norm_fc)) %>%
  ungroup(SampleDescriptions)

plottingTable2 = allData %>%
  dplyr::select(SampleDescriptions, Rcolonies_norm_fc, meanRcolonies, sdRcolonies, seRcolonies, meanRcellsOutsideColonies, sdRcellsOutsideColonies) %>%
  distinct()

plottingTable = allData %>%
  dplyr::select(SampleDescriptions, meanRcolonies, sdRcolonies, seRcolonies, meanRcellsOutsideColonies, sdRcellsOutsideColonies) %>%
  distinct()

plottingTable$SampleDescriptions = factor(plottingTable$SampleDescriptions,levels(plottingTable$SampleDescriptions)[c(2, 1, 4, 3)])

#PLot results
timingPlot_Rcolonies = ggplot(plottingTable, aes(x = reorder(SampleDescriptions, levels(plottingTable$SampleDescriptions)), y = meanRcolonies)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = (meanRcolonies - seRcolonies), ymax = (meanRcolonies + seRcolonies)), width = 0.2) +
  geom_point(data = plottingTable2, aes(x = SampleDescriptions, y = Rcolonies_norm_fc), stat = "identity", size = 2) +
  ylab("Number of resistant colonies (Fold change over control)") +
  xlab("Treatment regimen") +
  theme_classic()
timingPlot_Rcolonies
ggsave("04_plots/180817_ColonyGrowth31_EPZ5676_Timing/200831_timingExperiment_Rcolonies_rep1.PDF")

