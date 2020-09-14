#The goal of this script is to plot the colony growth data from of different KOs at different drug concentrations

#Read in data
baseline = read.csv("03_extractedData/200302_PLXdosing_limitedSet_responseCurve/baseline/colonyGrowthResults_PLX_DrugResponseCurve_baseline.csv", header = TRUE)
plateA_1uM = read.csv("03_extractedData/200302_PLXdosing_limitedSet_responseCurve/plateA_1uM/colonyGrowthResults_PLX_DrugResponseCurve_plateA_1uM.csv", header = TRUE)
plateA_2uM = read.csv("03_extractedData/200302_PLXdosing_limitedSet_responseCurve/plateA_2uM/colonyGrowthResults_PLX_DrugResponseCurve_plateA_2uM.csv", header = TRUE)
plateA_4uM = read.csv("03_extractedData/200302_PLXdosing_limitedSet_responseCurve/plateA_4uM/colonyGrowthResults_PLX_DrugResponseCurve_plateA_4uM.csv", header = TRUE)

plateB_1uM = read.csv("03_extractedData/200302_PLXdosing_limitedSet_responseCurve/plateB_1uM/colonyGrowthResults_PLX_DrugResponseCurve_plateB_1uM.csv", header = TRUE)
plateB_2uM = read.csv("03_extractedData/200302_PLXdosing_limitedSet_responseCurve/plateB_2uM/colonyGrowthResults_PLX_DrugResponseCurve_plateB_2uM.csv", header = TRUE)
plateB_4uM = read.csv("03_extractedData/200302_PLXdosing_limitedSet_responseCurve/plateB_4uM/colonyGrowthResults_PLX_DrugResponseCurve_plateB_4uM.csv", header = TRUE)

plateC_1uM = read.csv("03_extractedData/200302_PLXdosing_limitedSet_responseCurve/plateC_1uM/colonyGrowthResults_PLX_DrugResponseCurve_plateC_1uM.csv", header = TRUE)
plateC_2uM = read.csv("03_extractedData/200302_PLXdosing_limitedSet_responseCurve/plateC_2uM/colonyGrowthResults_PLX_DrugResponseCurve_plateC_2uM.csv", header = TRUE)
plateC_4uM = read.csv("03_extractedData/200302_PLXdosing_limitedSet_responseCurve/plateC_4uM/colonyGrowthResults_PLX_DrugResponseCurve_plateC_4uM.csv", header = TRUE)

#Add drug concentration
plateA_1uM$drugConcentration = "1uM"
plateA_2uM$drugConcentration = "2uM"
plateA_4uM$drugConcentration = "4uM"

plateB_1uM$drugConcentration = "1uM"
plateB_2uM$drugConcentration = "2uM"
plateB_4uM$drugConcentration = "4uM"

plateC_1uM$drugConcentration = "1uM"
plateC_2uM$drugConcentration = "2uM"
plateC_4uM$drugConcentration = "4uM"

#Add plate name
plateA_1uM$plateName = "plateA"
plateA_2uM$plateName = "plateA"
plateA_4uM$plateName = "plateA"

plateB_1uM$plateName = "plateB"
plateB_2uM$plateName = "plateB"
plateB_4uM$plateName = "plateB"

plateC_1uM$plateName = "plateC"
plateC_2uM$plateName = "plateC"
plateC_4uM$plateName = "plateC"

#Add sample name
baseline$sampleName = c("Neg01", "DOT1L_KO", "LATS2_KO", "BRD2_KO")

plateA_1uM$sampleName = c("Neg01", "DOT1L_KO", "LATS2_KO", "BRD2_KO")
plateA_2uM$sampleName = c("Neg01", "DOT1L_KO", "LATS2_KO", "BRD2_KO")
plateA_4uM$sampleName = c("Neg01", "DOT1L_KO", "LATS2_KO", "BRD2_KO")

plateB_1uM$sampleName = c("Neg01", "DOT1L_KO", "LATS2_KO", "BRD2_KO")
plateB_2uM$sampleName = c("Neg01", "DOT1L_KO", "LATS2_KO", "BRD2_KO")
plateB_4uM$sampleName = c("Neg01", "DOT1L_KO", "LATS2_KO", "BRD2_KO")

plateC_1uM$sampleName = c("Neg01", "DOT1L_KO", "LATS2_KO", "BRD2_KO")
plateC_2uM$sampleName = c("Neg01", "DOT1L_KO", "LATS2_KO", "BRD2_KO")
plateC_4uM$sampleName = c("Neg01", "DOT1L_KO", "LATS2_KO", "BRD2_KO")

#Joint data
jointData = bind_rows(plateA_1uM, plateA_2uM, plateA_4uM, plateB_1uM, plateB_2uM, plateB_4uM, plateC_1uM, plateC_2uM, plateC_4uM)

#Select useful fields
jointData = jointData %>%
  dplyr::select(well_number, sampleName, plateName, drugConcentration, num_colonies)

baseline = baseline %>%
  dplyr::select(sampleName, num_cells)

#Add baseline number of cells
jointData = left_join(jointData, baseline)

#Normalize number of colonies
jointData = jointData %>%
  mutate(num_colonies_norm = num_colonies * 10000 / num_cells)

#Obtain means and SE
jointData = jointData %>%
  group_by(sampleName, drugConcentration) %>%
  mutate(meanColonies = mean(num_colonies_norm)) %>%
  mutate(seColonies = sd(num_colonies_norm) / sqrt(3)) %>%
  ungroup(sampleName, drugConcentration)

#Make table for plotting
plottingTable = jointData %>%
  dplyr::select(well_number, sampleName, drugConcentration, meanColonies, seColonies) %>%
  distinct()

#Make plots: by drug concentration
coloniesByConcentration = ggplot(plottingTable, aes(x = drugConcentration, y = meanColonies, fill = as.character(well_number))) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = meanColonies - seColonies, ymax =  meanColonies + seColonies), width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = jointData, aes(x = drugConcentration, y = num_colonies_norm), stat = "identity", size = 1.5, position = position_dodge(0.9)) +
  geom_text(aes(x = drugConcentration, y = meanColonies, label = sampleName), position = position_dodge(0.9), angle = 90, size = 3) + 
  ylab("Number of resistant colonies (normalized to an input of 10,000 cells)") +
  xlab("PLX4032 concentration") +
  theme_classic()
coloniesByConcentration
ggsave("04_plots/200302_PLXdosing_limitedSet_responseCurve/200831_colonyGrowth_drugResponse_byConcentration.PDF")

# coloniesByConcentration = ggplot(plottingTable, aes(x = as.character(well_number), y = meanColonies, fill = drugConcentration)) +
#   geom_bar(stat = "identity", position = position_dodge(0.9)) +
#   geom_errorbar(aes(ymin = meanColonies - seColonies, ymax =  meanColonies + seColonies), width = 0.2, position = position_dodge(0.9)) +
#   geom_text(aes(x = as.character(well_number), y = meanColonies, label = sampleName), position = position_dodge(0.9), angle = 90, size = 3) + 
#   ylab("Number of resistant colonies (normalized to an input of 10,000 cells)") +
#   xlab("Knockout") +
#   theme_classic()
# coloniesByConcentration
# ggsave("04_plots/200302_PLXdosing_limitedSet_responseCurve/colonyGrowth_drugResponse_byKO.PDF")





