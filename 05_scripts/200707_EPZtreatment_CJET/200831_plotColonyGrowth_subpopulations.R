# The goal of this script is to plot the colony growth data from WM989 cells pre-treated with either DMSO or EPZ5676 (4uM) for seven days, then with
# PLX4032 1uM for 2-3 weeks.

#Read in data
baseline = read.csv("03_extractedData/200707_EPZtreatment_CJET/200707_colonyGrowthResults_EPZ5676treatment_SortedCells_baseline.csv", header = TRUE)
plateA = read.csv("03_extractedData/200707_EPZtreatment_CJET/200707_colonyGrowthResults_EPZ5676treatment_SortedCells_plateA.csv", header = TRUE)
plateB = read.csv("03_extractedData/200707_EPZtreatment_CJET/200707_colonyGrowthResults_EPZ5676treatment_SortedCells_plateB.csv", header = TRUE)
plateC = read.csv("03_extractedData/200707_EPZtreatment_CJET/200707_colonyGrowthResults_EPZ5676treatment_SortedCells_plateC.csv", header = TRUE)

#Add info to samples
baseline$plate = "baseline"
plateA$plate = "plateA"
plateB$plate = "plateB"
plateC$plate = "plateC"

baseline$sampleName = c("DMSO_bulk", "DMSO_low", "DMSO_high", "EPZ_bulk", "EPZ_low", "EPZ_high")
plateA$sampleName = c("DMSO_bulk", "DMSO_low", "DMSO_high", "EPZ_bulk", "EPZ_low", "EPZ_high")
plateB$sampleName = c("DMSO_bulk", "DMSO_low", "DMSO_high", "EPZ_bulk", "EPZ_low", "EPZ_high")
plateC$sampleName = c("DMSO_bulk", "DMSO_low", "DMSO_high", "EPZ_bulk", "EPZ_low", "EPZ_high")

baseline$preTreatment = c("DMSO", "DMSO", "DMSO", "EPZ", "EPZ", "EPZ")
plateA$preTreatment = c("DMSO", "DMSO", "DMSO", "EPZ", "EPZ", "EPZ")
plateB$preTreatment = c("DMSO", "DMSO", "DMSO", "EPZ", "EPZ", "EPZ")
plateC$preTreatment = c("DMSO", "DMSO", "DMSO", "EPZ", "EPZ", "EPZ")


baseline$phenotype = c("bulk", "low", "high", "bulk", "low", "high")
plateA$phenotype = c("bulk", "low", "high", "bulk", "low", "high")
plateB$phenotype = c("bulk", "low", "high", "bulk", "low", "high")
plateC$phenotype = c("bulk", "low", "high", "bulk", "low", "high")


#Keep Columns of interest
plateA = plateA %>%
  dplyr::select(plate, sampleName, preTreatment, phenotype, num_colonies)

plateB = plateB %>%
  dplyr::select(plate, sampleName, preTreatment, phenotype, num_colonies)

plateC = plateC %>%
  dplyr::select(plate, sampleName, preTreatment, phenotype, num_colonies)

baseline = baseline %>%
  dplyr::select(sampleName, num_cells)

#Join data into single table
jointData = bind_rows(plateA, plateB, plateC)
jointData = left_join(jointData, baseline, by = "sampleName")

#Obtain normalize number of colonies to 1,000 initial cells
jointData = jointData %>%
  mutate(num_colonies_norm = (num_colonies * 1000) / num_cells)

#Obtain mean and SE mean
jointData = jointData %>%
  group_by(sampleName) %>%
  mutate(mean_num_colonies_norm = mean(num_colonies_norm)) %>%
  mutate(sd_num_colonies_norm = sd(num_colonies_norm)) %>%
  mutate(se_num_colonies_norm = sd(num_colonies_norm)/sqrt(3)) %>%
  ungroup(sampleName)

#Make plotting table
plottingTable = jointData %>%
  dplyr::select(sampleName, preTreatment, phenotype, num_colonies_norm,mean_num_colonies_norm, sd_num_colonies_norm, se_num_colonies_norm)
  
#Plot means and SD, with dot plot
graphOfColonies = ggplot(plottingTable, aes(x = phenotype, y = mean_num_colonies_norm, fill = preTreatment)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = (mean_num_colonies_norm - sd_num_colonies_norm), ymax = (mean_num_colonies_norm + sd_num_colonies_norm)), width = 0.2, position = position_dodge(0.9)) +
  ylab("Number of resistant colonies (normalized to an input of 1,000 cells)") +
  geom_point(data = plottingTable, aes(x = phenotype, y = num_colonies_norm), stat = "identity", position = position_dodge(0.9), size = 2) +
  xlab("Phenotype") +
  theme_classic()
graphOfColonies
ggsave("04_plots/200707_EPZtreatment_CJET/200707_colonyGrowth_EPZ5676vsDMSO_subpopulations.PDF")

















