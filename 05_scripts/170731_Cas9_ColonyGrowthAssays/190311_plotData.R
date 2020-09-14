#The goal of this script is to plot colony growth data from WM989 and WM989-cas9 NGFR-EGFR sorts
#Read in summarized results
wt_DP = read.csv("01_rawData/170731_Cas9_ColonyGrowthAssays/170719_WM989A6G3_DP/colonyGrowthResults_WM989_DP.csv", header = T)
wt_NGFR = read.csv("01_rawData/170731_Cas9_ColonyGrowthAssays/170719_WM989A6G3_NGFR/colonyGrowthResults_WM989_NGFR.csv", header = T)
wt_EGFR = read.csv("01_rawData/170731_Cas9_ColonyGrowthAssays/170719_WM989A6G3_EGFR/colonyGrowthResults_WM989_EGFR.csv", header = T)
wt_ungated = read.csv("01_rawData/170731_Cas9_ColonyGrowthAssays/170719_WM989A6G3_mix/colonyGrowthResults_WM989_ungated.csv", header = T)

cas9_DP = read.csv("01_rawData/170731_Cas9_ColonyGrowthAssays/170721_WM989cas95a3_DP/colonyGrowthResults_WM989cas9_DP.csv", header = T)
cas9_NGFR = read.csv("01_rawData/170731_Cas9_ColonyGrowthAssays/170722_WM989cas95a3_NGFR/colonyGrowthResults_WM989cas9_NGFR.csv", header = T)
cas9_EGFR = read.csv("01_rawData/170731_Cas9_ColonyGrowthAssays/170722_WM989cas95a3_EGFR/colonyGrowthResults_WM989cas9_EGFR.csv", header = T)
cas9_ungated = read.csv("01_rawData/170731_Cas9_ColonyGrowthAssays/170721_WM989cas95a3_mix/colonyGrowthResults_WM989cas9_ungated.csv", header = T)

#Add sample info
wt_DP = wt_DP %>%
  mutate(sampleName = "DP") %>%
  mutate(cellLines = "WM989")
wt_NGFR = wt_NGFR %>%
  mutate(sampleName = "NGFR") %>%
  mutate(cellLines = "WM989")
wt_EGFR = wt_EGFR %>%
  mutate(sampleName = "EGFR") %>%
  mutate(cellLines = "WM989")
wt_ungated = wt_ungated %>%
  mutate(sampleName = "ungated") %>%
  mutate(cellLines = "WM989")

cas9_DP = cas9_DP %>%
  mutate(sampleName = "DP") %>%
  mutate(cellLines = "WM989cas9")
cas9_NGFR = cas9_NGFR %>%
  mutate(sampleName = "NGFR") %>%
  mutate(cellLines = "WM989cas9")
cas9_EGFR = cas9_EGFR %>%
  mutate(sampleName = "EGFR") %>%
  mutate(cellLines = "WM989cas9")
cas9_ungated = cas9_ungated %>%
  mutate(sampleName = "ungated") %>%
  mutate(cellLines = "WM989cas9")

#Join all data
allData = bind_rows(wt_DP, wt_NGFR, wt_EGFR, wt_ungated, cas9_DP, cas9_NGFR, cas9_EGFR, cas9_ungated)

#Convert the sample name into factors
allData$sampleName = factor(allData$sampleName, levels = c("ungated", "EGFR", "NGFR", "DP"))

#Plot data
barplot_colonyGrowth = ggplot(data = allData, aes(y = num_colonies, x = sampleName, fill = sampleName)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(. ~ cellLines) +
  theme_classic()
barplot_colonyGrowth
ggsave("04_plots/170731_Cas9_ColonyGrowthAssays/barplot_WM989.PDF")
