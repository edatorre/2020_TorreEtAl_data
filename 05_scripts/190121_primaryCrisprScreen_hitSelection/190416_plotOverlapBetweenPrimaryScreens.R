#The goal of this script is to plot the results of the fate screen against the results from the state screen

#Load the screen results
screenResults = read.table("02_metadata/190402_metadataOfHits.txt", header = T)

#Separate hits from the rest of the population
DPhits = screenResults %>%
  filter(DPtier == "tier075" | DPtier == "tier066") %>%
  dplyr::select(geneName)

VemRhits = screenResults %>%
  filter(VemRtier == "tier075" | VemRtier == "tier066") %>%
  dplyr::select(geneName)

doubleHits = DPhits %>%
  filter(geneName %in% VemRhits$geneName) %>%
  dplyr::select(geneName)

#Make list of all hits and add group labels
hits = bind_rows(DPhits, VemRhits)
hits = hits %>%
  distinct()

hits$group = "none"
hits$group = ifelse(hits$geneName %in% DPhits$geneName, "DP", hits$group)
hits$group = ifelse(hits$geneName %in% VemRhits$geneName, "VemR", hits$group)
hits$group = ifelse(hits$geneName %in% doubleHits$geneName, "both", hits$group)


#Make plottingTable of hits
hits = left_join(hits, screenResults, by = "geneName")

#Make table for names
doubleHits = left_join(doubleHits, screenResults, by = "geneName")


#Plot scatterplot
plotScreenResults = ggplot() +
  geom_point(data = hits, aes(x = hits$medianlFC_DP, y = hits$medianlFC_VemR, color = hits$group)) +
  geom_density_2d(data = screenResults, aes(x = screenResults$medianlFC_DP, y = screenResults$medianlFC_VemR)) +
  stat_density_2d(data = screenResults, aes(x = screenResults$medianlFC_DP, y = screenResults$medianlFC_VemR, fill = ..level..), geom = "polygon") +
  geom_text_repel(data = doubleHits, aes(x = doubleHits$medianlFC_DP, y = doubleHits$medianlFC_VemR, label = geneName)) +
  scale_color_manual(values = c("both" = "purple", "DP" = "orange", "VemR" = "blue")) +
  scale_fill_gradientn(colours = brewer.pal( 7, "Greys" )) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_classic()
plotScreenResults
ggsave("04_plots/190121_primaryCrisprScreen_hitSelection/overlapBetweenScreens.PDF")
