#The goal of this script is to plot the negative seleciton assay
#Input data
allData = read.table("01_rawData/170215_WM989cas9_NegativeSelectionAssay/negativeSelection_3_allData.txt", header = TRUE)

#make tall table
allData_tall = allData %>%
  gather(day, percentPositive, 4:10)

#Eliminate null control form data
allData_tall = allData_tall %>%
  filter(gRNA != "Null")

#obtain normalization values
normalizationValues = allData_tall %>%
  filter(day == "Day5")

normalizationValues = normalizationValues %>%
  group_by(CellLine, gRNA) %>%
  mutate(normalizationFactor = mean(percentPositive)) %>%
  ungroup(CellLine, gRNA) %>%
  dplyr::select(CellLine, gRNA, normalizationFactor) %>%
  distinct()

#Add normalization values to allData
allData_tall = left_join(allData_tall, normalizationValues, by = c("CellLine", "gRNA"))

#Normalize data to day 5
allData_tall = allData_tall %>%
  mutate(percentPositive_norm = percentPositive * 100 / normalizationFactor)

#Now obtain means and standard deviations
allData_tall = allData_tall %>%
  group_by(CellLine, gRNA, day) %>%
  mutate(meanPercentPositive = mean(percentPositive_norm)) %>%
  mutate(sdPercentPositive = sd(percentPositive_norm)) %>%
  ungroup(CellLine, gRNA, day)

allData_summarized = allData_tall %>%
  dplyr::select(CellLine, gRNA, day, meanPercentPositive, sdPercentPositive) %>%
  distinct()

allData_summarized[is.na(allData_summarized)] = 0

allData_summarized$day = factor(allData_summarized$day)
allData_summarized$gRNA = factor(allData_summarized$gRNA)

allData_summarized$gRNA = factor(allData_summarized$gRNA, levels = c("ROSA", "CDK9", "RPA2", "Null"))

allData_summarized$day = factor(allData_summarized$day, levels = c("Day5", "Day9", "Day14", "Day21", "Day28", "Day36", "Day43"))


#Plot data
negativeSelectionPlot = ggplot(data = allData_summarized, aes(x = gRNA, y = meanPercentPositive, fill = day)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(CellLine ~ .) +
  theme_classic() +
  ylab("Percent of total population (%GFP+)")
negativeSelectionPlot
ggsave("04_plots/170215_WM989cas9_NegativeSelectionAssay/negativeSelection_barplots.PDF")



