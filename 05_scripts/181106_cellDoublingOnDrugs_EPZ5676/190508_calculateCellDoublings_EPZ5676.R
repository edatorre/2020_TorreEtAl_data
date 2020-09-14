# #The goal of this script is to process cell doubling data

#Read input data
day1 = read.table("01_rawData/181106_cellDoublingOnDrugs_EPZ5676/day1.txt")
day3 = read.table("01_rawData/181106_cellDoublingOnDrugs_EPZ5676/day3.txt")
day5 = read.table("01_rawData/181106_cellDoublingOnDrugs_EPZ5676/day5.txt")
day7 = read.table("01_rawData/181106_cellDoublingOnDrugs_EPZ5676/day7.txt")

#Load plate map
plateMap = read.table("02_metadata/181106_cellDoublingOnDrugs_EPZ5676/plateMap.txt", header = T)

#Assign column and row names to data
columnNames = c("rep1", "rep2", "rep3", "background")
rowNames = c("row1", "row2", "row3")

colnames(day1) = columnNames
colnames(day3) = columnNames
colnames(day5) = columnNames
colnames(day7) = columnNames

rownames(day1) = rowNames
rownames(day3) = rowNames
rownames(day5) = rowNames
rownames(day7) = rowNames

#Make tall tables
plateMap$row = row.names(plateMap)
day1$row = row.names(day1)
day3$row = row.names(day3)
day5$row = row.names(day5)
day7$row = row.names(day7)

plateMap = plateMap %>%
  gather("rep", "sample", 1:3)

day1 = day1 %>%
  gather("rep", "value", 1:3)

day3 = day3 %>%
  gather("rep", "value", 1:3)

day5 = day5 %>%
  gather("rep", "value", 1:3)

day7 = day7 %>%
  gather("rep", "value", 1:3)

#Subtract backgorund and add smaple names
day1 = day1 %>%
  mutate(d1_value_norm = value - background) %>%
  left_join(., plateMap, by = c("row", "rep")) %>%
  dplyr::select(sample, rep, d1_value_norm)

day3 = day3 %>%
  mutate(d3_value_norm = value - background) %>%
  left_join(., plateMap, by = c("row", "rep")) %>%
  dplyr::select(sample, rep, d3_value_norm)

day5 = day5 %>%
  mutate(d5_value_norm = value - background) %>%
  left_join(., plateMap, by = c("row", "rep")) %>%
  dplyr::select(sample, rep, d5_value_norm)

day7 = day7 %>%
  mutate(d7_value_norm = value - background) %>%
  left_join(., plateMap, by = c("row", "rep")) %>%
  dplyr::select(sample, rep, d7_value_norm)

#Join data
allData = left_join(day1, day3, by = c("sample", "rep")) %>%
  left_join(., day5, by = c("sample", "rep")) %>%
  left_join(., day7, by = c("sample", "rep"))

#Add number of samples 
allData = allData %>%
  group_by(sample) %>%
  mutate(nSamples = n_distinct(rep)) %>%
  ungroup(sample)

#limit dataset to EPZ5676
allData = allData %>%
  filter(sample != "PLX4032_1uM")

#Make plotting table
plottingTable = allData %>%
  dplyr::select(sample, rep, nSamples, d1_value_norm, d3_value_norm, d5_value_norm, d7_value_norm) %>%
  gather("timepoint", "measurement", 4:7)

plottingTable = plottingTable %>%
  group_by(sample, timepoint) %>%
  mutate(meanMeasurement = mean(measurement)) %>%
  mutate(sdMeasurement = sd(measurement)) %>%
  mutate(seMeasurement = sdMeasurement / sqrt(nSamples))

plottingTable = plottingTable %>%
  dplyr::select(-rep, -measurement, -nSamples) %>%
  distinct()

growthPlot = ggplot(data = plottingTable, aes(x = timepoint, y = meanMeasurement, group = sample, color = sample)) +
  geom_errorbar(aes(ymin = meanMeasurement - seMeasurement, ymax = meanMeasurement + seMeasurement), width=.1) +
  geom_line() +
  geom_point() +
  theme_classic()
growthPlot
ggsave("04_plots/181106_cellDoublingOnDrugs_EPZ5676/cellDoubling_EPZ5676_lineGraph.PDF")
