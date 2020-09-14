# #The goal of this script is to create a single melted couts file

#load files
dataset1 = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/A6_ColonyE3_NoDrug/htseq/A6_ColonyE3_NoDrug.htseq.stdout", sep = "\t")
dataset2 = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/A6_ColonyG3_NoDrug/htseq/A6_ColonyG3_NoDrug.htseq.stdout", sep = "\t")
dataset3 = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/A6_ColonyH1_NoDrug/htseq/A6_ColonyH1_NoDrug.htseq.stdout", sep = "\t")
dataset4 = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/ColonyE3-WellC6/htseq/ColonyE3-WellC6.htseq.stdout", sep = "\t")
dataset5 = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/ColonyG3-WellB6-drug/htseq/ColonyG3-WellB6-drug.htseq.stdout", sep = "\t")
dataset6 = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/ColonyH1-WellC2/htseq/ColonyH1-WellC2.htseq.stdout", sep = "\t")
dataset7 = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/EGFRsorts/meltedData_20150722_EGFRsortRNArun1.tsv", sep = "\t", header = T)
dataset8 = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/EGFRsorts/meltedData_20150831_EGFRsortRNArun2.tsv", sep = "\t", header = T)
dataset9 = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/ColonyG3-WellC6-drug/htseq/ColonyG3-WellC6-drug.htseq.stdout", sep = "\t")
dataset10 = read.table("01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/NGFRsort_ben/meltedData.tsv", sep = "\t", header = T)

#Add sample name to each file
dataset1 = dataset1 %>%
  mutate(sampleID = "A6_ColonyE3_NoDrug")

dataset2 = dataset2 %>%
  mutate(sampleID = "A6_ColonyG3_NoDrug")

dataset3 = dataset3 %>%
  mutate(sampleID = "A6_ColonyH1_NoDrug")

dataset4 = dataset4 %>%
  mutate(sampleID = "ColonyE3-WellC6")

dataset5 = dataset5 %>%
  mutate(sampleID = "ColonyG3-WellB6-drug")

dataset6 = dataset6 %>%
  mutate(sampleID = "ColonyH1-WellC2")

dataset9 = dataset9 %>%
  mutate(sampleID = "ColonyG3-WellC6-drug")

#merge data
allData_1 = bind_rows(dataset1, dataset2, dataset3, dataset4, dataset5, dataset6, dataset9)

#rename columns
colnames(allData_1) = c("gene_ID", "counts", "sampleID")

#merge data
allData_2 = bind_rows(dataset7, dataset8)
allData_2$experiment = NULL
colnames(allData_2) = c("sampleID", "gene_ID", "counts")

#Work on Ben's data
dataset10$experiment = NULL
colnames(dataset10) = c("sampleID", "gene_ID", "counts")

#Merge all data
allData = bind_rows(allData_1, allData_2, dataset10)

#save table
write.table(allData, file = "01_rawData/190408_RNAseq_otherSamples/otherDatasets/rawCounts/allCounts_meltedData.txt", sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)

