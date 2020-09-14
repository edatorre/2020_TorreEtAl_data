# #make plot of sliding scale

# #set working directory
# setwd("/Users/etorre/Dropbox (RajLab)/Eduardo_shared/et_paperCRISPR/data/")

#Load name conversion file
nameKey = read.table("02_metadata/190402_screenTargets_geneDomainTarget_nameConversion.txt", header = T)

#Define the positive controls
positiveControls = c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
negativeControls = c("neg")

#Load IF validation
CGvalidation = read.table("03_extractedData/190406_validationTables/validationTable_VemR.txt", header = T)
CGvalidation = CGvalidation %>%
  dplyr::select(geneName, CG_validationStatus)

#Load screen results and merge them
epigenetic_screenResults_VemR = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_epigenetic_allInfo.txt", header = T, sep = "\t")
kinase_screenResults_VemR = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_kinase_allInfo.txt", header = T, sep = "\t")
tf_screenResults_VemR = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemRanalysis_TF_allInfo.txt", header = T, sep = "\t")

screenResults_VemR = bind_rows(epigenetic_screenResults_VemR, kinase_screenResults_VemR, tf_screenResults_VemR)

#Eliminate controls from these data
screenResults_VemR = screenResults_VemR %>%
  dplyr::select(gene, domain, lFC, medianlFC) %>%
  filter(!gene %in% positiveControls) %>%
  filter(!gene %in% negativeControls)

#Load sliding scale tables
epigenetic_slidingScaleTable = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemR_epigeneticLibrary_slidingScaleTable.txt", header = T)
kinase_slidingScaleTable = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemR_kinaseLibrary_slidingScaleTable.txt", header = T)
tf_slidingScaleTable = read.table("03_extractedData/190121_primaryCrisprScreen_hitSelection/VemR_TFLibrary_slidingScaleTable.txt", header = T)

#combine sliding scale tables
slidingScaleTable = bind_rows(epigenetic_slidingScaleTable, kinase_slidingScaleTable, tf_slidingScaleTable)
slidingScaleTable = slidingScaleTable %>%
  dplyr::select(-gRNA) %>%
  distinct()

#Combine the screen restuls with the sliding scale table
screenResults_VemR = inner_join(screenResults_VemR, slidingScaleTable, by = "domain")
screenResults_VemR = screenResults_VemR %>%
  distinct()

#Add in validation data
CGvalidation$geneName = gsub("BRD2-BD1", "BRD2", CGvalidation$geneName) #Here I fix the BRD2 name to be able to combine the tables
screenResults_VemR = left_join(screenResults_VemR, CGvalidation, by = c("gene" = "geneName")) #Now combine the data

#break tables into sections
section1 = screenResults_VemR %>% 
  distinct() %>%
  filter(hitStatus_075 == "yes") 
breakpoint_secrtion1 = nrow(section1)
section1 = section1 %>%
  mutate(section = "section1_75percent") %>%
  mutate(filter = "75%") %>%
  mutate(plottingID = 1:breakpoint_secrtion1)

section2 = screenResults_VemR %>%
  distinct() %>%
  filter(hitStatus_066 == "yes") %>%
  filter(hitStatus_075 == "No")
breakpoint_secrtion2 = (breakpoint_secrtion1 + 19) + nrow(section2)
section2 = section2 %>%
  mutate(section = "section2_66percent") %>%
  mutate(filter = "66%") %>%
  mutate(plottingID = (breakpoint_secrtion1 + 20):(breakpoint_secrtion2))

section3 = screenResults_VemR %>%
  distinct() %>%
  filter(hitStatus_050 == "yes") %>%
  filter(hitStatus_066 == "No")
breakpoint_secrtion3 = (breakpoint_secrtion2 + 19) + nrow(section3)
section3 = section3 %>%
  mutate(section = "section3_50percent") %>%
  mutate(filter = "50%") %>%
  mutate(plottingID = (breakpoint_secrtion2 + 20):(breakpoint_secrtion3))

section4 = screenResults_VemR %>%
  #dplyr::select(-gRNA) %>%
  distinct() %>%
  filter(hitStatus_050 == "No")
breakpoint_secrtion4 = (breakpoint_secrtion3 + 19) + nrow(section4)
section4 = section4 %>%
  mutate(section = "section4_under50percent") %>%
  mutate(filter = "none") %>%
  mutate(plottingID = (breakpoint_secrtion3 + 20):(breakpoint_secrtion4))

#Merge all sections and all breakpoints
plottingTable = bind_rows(section1, section2, section3, section4)

#Remove controls from the table (in case one of the previous tables added them back)
plottingTable = plottingTable %>%
  filter(!domain %in% negativeControls) %>%
  filter(!domain %in% positiveControls)

#calculate number of domains per "tier"
plottingTable = plottingTable %>%
  group_by(section) %>%
  mutate(nDomain = n_distinct(domain)) %>%
  ungroup(section)

#calculate the number of domans tested per domain
tmp = plottingTable %>%
  filter(CG_validationStatus == "tested_validated" | CG_validationStatus == "tested_NotValidated") %>%
  dplyr::select(domain, section, CG_validationStatus) %>%
  group_by(section) %>%
  mutate(nHitsTested = n_distinct(domain)) %>%
  ungroup(section) %>%
  dplyr::select(section, nHitsTested) %>%
  distinct()

plottingTable = left_join(plottingTable, tmp, by = "section")

#Calulate the number of hits tested AND validated
tmp = plottingTable %>%
  filter(CG_validationStatus == "tested_validated") %>%
  dplyr::select(domain, section, CG_validationStatus) %>%
  group_by(section) %>%
  mutate(nHitsValidated = n_distinct(domain)) %>%
  ungroup(section) %>%
  dplyr::select(section, nHitsValidated) %>%
  distinct()

plottingTable = left_join(plottingTable, tmp, by = "section")

#Calculate percent of hita validated (out of the ones tested)
plottingTable$nHitsTested[is.na(plottingTable$nHitsTested)] = 0 
plottingTable$nHitsValidated[is.na(plottingTable$nHitsValidated)] = 0  
plottingTable = plottingTable %>%
  mutate(validationRate = (nHitsValidated / nHitsTested) * 100)

#Make a table that contains the validation statistics for the graph
plottingTable_stats = plottingTable %>%
  dplyr::select(section, plottingID, filter, validationRate) %>%
  group_by(section) %>%
  mutate(plottingID_text_x = min(plottingID) + ((max(plottingID) - min(plottingID))/ 2)) %>%
  mutate(plottingID_text_y1 = -6) %>%
  mutate(plottingID_text_y2 = 6) %>%
  dplyr::select(-plottingID) %>%
  distinct()

#Make the plot
plot_slidingScale <- ggplot() +
  geom_point(data = plottingTable, shape = 20, aes(x = -plottingID, y = lFC, color = CG_validationStatus)) +
  scale_color_manual(values = c("notTested" = "gray", "tested_validated" = "red", "tested_NotValidated" = "black")) + 
  scale_size_manual(values = c("notTested" = 0.5, "tested_validated" = 2, "tested_NotValidated" = 2)) +
  geom_text(data = plottingTable_stats, aes(x = -plottingID_text_x, y = plottingID_text_y1, label = filter), size = 2) +
  geom_text(data = plottingTable_stats, aes(x = -plottingID_text_x, y = plottingID_text_y2, label = validationRate), size = 2) +
  theme(legend.position="none") +
  scale_y_continuous(breaks = round(seq(-20, 10, by = 1),1)) +
  geom_hline(yintercept = c(1, -1), linetype="dashed") +
  geom_vline(xintercept = c((-1 * breakpoint_secrtion1), (-1 * breakpoint_secrtion2), (-1 * breakpoint_secrtion3), (-1 * breakpoint_secrtion4))) +
  coord_flip() +
  theme_classic()
plot(plot_slidingScale)
ggsave("04_plots/190406_validationTables/VemR_validationPlot_byGuides.PDF", width = 5, height = 7)
