#The goal of this script is to plot the colony growth data from WM989 and WM989-cas9 cells

#Read in data
baseline = read.csv("03_extractedData/200302_PLXdosing_limitedSet_WM989comp/baseline/colonyGrowthResults_WM989_baseline.csv", header = TRUE)
plateA = read.csv("03_extractedData/200302_PLXdosing_limitedSet_WM989comp/plateA_1uM/colonyGrowthResults_WM989_plateA_1uM.csv", header = TRUE)
plateB = read.csv("03_extractedData/200302_PLXdosing_limitedSet_WM989comp/plateB_1uM/colonyGrowthResults_WM989_plateB_1uM.csv", header = TRUE)
plateC = read.csv("03_extractedData/200302_PLXdosing_limitedSet_WM989comp/plateC_1uM/colonyGrowthResults_WM989_plateC_1uM.csv", header = TRUE)

#Add info to samples
baseline$plate = "baseline"
plateA$plate = "plateA"
plateB$plate = "plateB"
plateC$plate = "plateC"

#Add sample name
baseline$sampleName = c("WM989", "WM989-cas9")
plateA$sampleName = c("WM989", "WM989-cas9")
plateB$sampleName = c("WM989", "WM989-cas9")
plateC$sampleName = c("WM989", "WM989-cas9")

#Keep Columns of interest
plateA = plateA %>%
  dplyr::select(plate, sampleName, num_colonies)

plateB = plateB %>%
  dplyr::select(plate, sampleName, num_colonies)

plateC = plateC %>%
  dplyr::select(plate, sampleName, num_colonies)

baseline = baseline %>%
  dplyr::select(sampleName, num_cells)

#Join data into single table
jointData = bind_rows(plateA, plateB, plateC)
jointData = left_join(jointData, baseline, by = "sampleName")

#Normalize number of colonies to 10,000 initial cells
jointData = jointData %>%
  mutate(num_colonies_norm = (num_colonies * 10000) / num_cells)

#Obtain mean
jointData = jointData %>%
  group_by(sampleName) %>%
  mutate(mean_colonies_norm = mean(num_colonies_norm)) %>%
  ungroup(sampleName)

#Obtain SE mean
jointData = jointData %>%
  group_by(sampleName) %>%
  mutate(se_colonies_norm = sd(num_colonies_norm) / sqrt(3)) %>% # samples in triplicates. n = 3
  ungroup(sampleName)

#make ploting table
plottingTable = jointData %>%
  dplyr::select(sampleName, mean_colonies_norm, se_colonies_norm) %>%
  distinct()

#plot the results
graphOfColonies = ggplot(plottingTable, aes(x = sampleName, y = mean_colonies_norm)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = (mean_colonies_norm - se_colonies_norm), ymax = (mean_colonies_norm + se_colonies_norm)), width = 0.2) +
  geom_point(data = jointData, aes(x = sampleName, y = num_colonies_norm), stat = "identity", size = 2, position = position_dodge(0.9)) +
  ylab("Number of resistant colonies (normalized to an input of 10,000 cells)") +
  xlab("Cell line") +
  theme_classic()
graphOfColonies
ggsave("04_plots/200302_PLXdosing_limitedSet_WM989comp/200831_colonyGrowth_WM989comparison.PDF")
