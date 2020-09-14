#The goal of this script is to produce a table with the number of mice in each group and the pValue for each comparison

#Load reformated data
allData = read.table("03_extractedData/190824_inVivoAnalysis/reformatedData.txt", header = T, sep = "\t")

#Count the number of mice
allData = allData %>%
  group_by(KO, treatmentArm, measurementGroup) %>%
  mutate(nMice = n_distinct(mouseID)) %>%
  ungroup(KO, treatmentArm, measurementGroup)

#Make summarize data table
summarizedData = allData %>%
  dplyr::select(KO, treatmentArm, measurementGroup, nMice) %>%
  distinct()

#Add number of control mice
KO_summarizedData = summarizedData %>%
  filter(KO != "Neg")

Control_summarizedData = summarizedData %>%
  filter(KO == "Neg")

#Change column names
colnames(KO_summarizedData)[4] = "nMice_KO"
colnames(Control_summarizedData)[4] = "nMice_Control"

#Keep only relevant columns
Control_summarizedData = Control_summarizedData %>%
  dplyr::select(-KO)

#join Data
KO_summarizedData = left_join(KO_summarizedData, Control_summarizedData, by = c("treatmentArm", "measurementGroup"))

#Load p values - BRD2
BRD2_treated = read.table("03_extractedData/190824_inVivoAnalysis/ttests_BRD2_treated_less.txt", header = T, sep = "\t")
BRD2_treated$measurementGroup = rownames(BRD2_treated)
BRD2_treated$KO = "BRD2"
BRD2_treated$treatmentArm = "treated"
colnames(BRD2_treated)[1] = "pVal"

BRD2_untreated = read.table("03_extractedData/190824_inVivoAnalysis/ttests_BRD2_untreated_less.txt", header = T, sep = "\t")
BRD2_untreated$measurementGroup = rownames(BRD2_untreated)
BRD2_untreated$KO = "BRD2"
BRD2_untreated$treatmentArm = "untreated"
colnames(BRD2_untreated)[1] = "pVal"

#Load p values - DOT1L
DOT1L_treated = read.table("03_extractedData/190824_inVivoAnalysis/ttests_DOT1L_treated_greater.txt", header = T, sep = "\t")
DOT1L_treated$measurementGroup = rownames(DOT1L_treated)
DOT1L_treated$KO = "DOT1L"
DOT1L_treated$treatmentArm = "treated"
colnames(DOT1L_treated)[1] = "pVal"

DOT1L_untreated = read.table("03_extractedData/190824_inVivoAnalysis/ttests_DOT1L_untreated_greater.txt", header = T, sep = "\t")
DOT1L_untreated$measurementGroup = rownames(DOT1L_untreated)
DOT1L_untreated$KO = "DOT1L"
DOT1L_untreated$treatmentArm = "untreated"
colnames(DOT1L_untreated)[1] = "pVal"

#Load p values - LATS2
LATS2_treated = read.table("03_extractedData/190824_inVivoAnalysis/ttests_LATS2_treated_greater.txt", header = T, sep = "\t")
LATS2_treated$measurementGroup = rownames(LATS2_treated)
LATS2_treated$KO = "LATS2"
LATS2_treated$treatmentArm = "treated"
colnames(LATS2_treated)[1] = "pVal"

LATS2_untreated = read.table("03_extractedData/190824_inVivoAnalysis/ttests_LATS2_untreated_greater.txt", header = T, sep = "\t")
LATS2_untreated$measurementGroup = rownames(LATS2_untreated)
LATS2_untreated$KO = "LATS2"
LATS2_untreated$treatmentArm = "untreated"
colnames(LATS2_untreated)[1] = "pVal"

#Combine P values
pValues = bind_rows(BRD2_treated, BRD2_untreated, DOT1L_treated, DOT1L_untreated, LATS2_treated, LATS2_untreated)

#Combine p value with nMive data
combinedData = left_join(KO_summarizedData, pValues, by = c("KO", "treatmentArm", "measurementGroup"))

#Save table
write.table(combinedData, file = "03_extractedData/190824_inVivoAnalysis/pValuesAndNmice.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
