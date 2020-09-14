# #Plot all growth curves.

#Read data
allData = read.table("03_extractedData/190824_inVivoAnalysis/reformatedData.txt", header = TRUE, sep = "\t")

medianTumorVolume = allData %>%
  group_by(treatmentArm, KO, measurementGroup) %>%
  mutate(meadianTumorVolume = median(tumorVolume)) %>%
  ungroup(treatmentArm, KO, measurementGroup) %>%
  dplyr::select(-mouseID, -tumorVolume) %>%
  distinct()


#BRD2 plot
BRD2_plottingTable = allData %>%
  filter(KO == "Neg" | KO == "BRD2")

BRD2_growthCurves = ggplot() +
  geom_point(data = BRD2_plottingTable, aes(x = treatmentDuration, y = tumorVolume, group = mouseID, color = KO)) +
  geom_line(data = BRD2_plottingTable, aes(x = treatmentDuration, y = tumorVolume, group = mouseID, color = KO), alpha = 0.4) +
  facet_grid(.~treatmentArm) + 
  theme_classic()
BRD2_growthCurves
ggsave("04_plots/190824_inVivoAnalysis/BRD2_allDataPoints.PDF")

#LATS2 plot
LATS2_plottingTable = allData %>%
  filter(KO == "Neg" | KO == "LATS2")

LATS2_growthCurves = ggplot() +
  geom_point(data = LATS2_plottingTable, aes(x = treatmentDuration, y = tumorVolume, group = mouseID, color = KO)) +
  geom_line(data = LATS2_plottingTable, aes(x = treatmentDuration, y = tumorVolume, group = mouseID, color = KO), alpha = 0.4) +
  facet_grid(.~treatmentArm) + 
  theme_classic()
LATS2_growthCurves
ggsave("04_plots/190824_inVivoAnalysis/LATS2_allDataPoints.PDF")

#DOT1L plot
DOT1L_plottingTable = allData %>%
  filter(KO == "Neg" | KO == "DOT1L")

DOT1L_growthCurves = ggplot() +
  geom_point(data = DOT1L_plottingTable, aes(x = treatmentDuration, y = tumorVolume, group = mouseID, color = KO)) +
  geom_line(data = DOT1L_plottingTable, aes(x = treatmentDuration, y = tumorVolume, group = mouseID, color = KO), alpha = 0.4) +
  facet_grid(.~treatmentArm) + 
  theme_classic()
DOT1L_growthCurves
ggsave("04_plots/190824_inVivoAnalysis/DOT1L_allDataPoints.PDF")
