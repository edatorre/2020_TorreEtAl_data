
runDESeq = function(geneOfInterest) {
  
  #select count data for analysis: controls
  t02_ControlSample = allCounts_rawCountsMatrix %>%
    select(one_of(c("neg01", "neg02", "neg03", "neg04", "neg05", "neg06", "neg07", "neg08", "neg09", "neg10")))
  
  #Add a column of row names to merge data
  t02_ControlSample = t02_ControlSample %>%
    mutate(gene_id = row.names(t02_ControlSample))
  
  # #I can use this line as a control to run the function within this script
  # geneOfInterest = "SOX10"
  
  #select count data for analysis: KOs
  t02_experimentalSample = allCounts_rawCountsMatrix %>%
    select(contains(geneOfInterest))
    
  #Add a column of row names to merge data
  t02_experimentalSample = t02_experimentalSample %>%
    mutate(gene_id = row.names(t02_experimentalSample))
  
  #Merge datasets
  t02_mergedData = left_join(t02_ControlSample, t02_experimentalSample, by = "gene_id")
  row.names(t02_mergedData) = t02_mergedData$gene_id
  t02_mergedData$gene_id = NULL

  #Here I remove empty columns or columns with very few reads
  readsPerSample = colSums(t02_mergedData)
  columnToKeep = which(readsPerSample > 500000) #Remove samples with less than 500K reads
  t02_mergedData = t02_mergedData[, c(columnToKeep)]
  
  #Make conditions list for DESeq. Has to be here because the number of columns changes later on.
  condition = c(rep("KO", ncol(t02_mergedData)))
  condition[1:10] = "wt"
  condition = factor(condition, levels = c("wt", "KO"))
  
  #Generate column data
  coldata = data.frame(condition)
  coldata$type = factor("paired-end")
  rownames(coldata) = colnames(t02_mergedData)
  
  #If there is at least one sample that is not a control then run loop. 
  if(ncol(t02_mergedData) > 10) {
    
    #Generate DEseq data set.
    dds <- DESeqDataSetFromMatrix(countData = t02_mergedData,
                                  colData = coldata,
                                  design = ~ condition)
    #Run DEseq
    dds <- DESeq(dds)
    
    #Extract Results
    res <- results(dds)
    res = data.frame(res)
    res$id = rownames(res)
    
    #Add KO info to results
    res = res %>%
      mutate(sampleKO = geneOfInterest)
    
    #filter to keep differentially expressed genes
    res_filtered = res %>%
      filter(padj < 0.05) %>%
      filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>%
      select(id, log2FoldChange, sampleKO)
    
    
    #save DESeq output
    write.table( res, file = paste("03_extractedData/190321_KO-RNAseq/DEseq/", geneOfInterest, "_differentialExpression_DESeq.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
    #save DESeq output in a single combined file
    write.table( res, file = "03_extractedData/190321_KO-RNAseq/DEseq/differentialExpression_DESeq_allTargets.txt", sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = TRUE)
    #Append the list of diffenretially expressed genes to an output file.
    write.table(res_filtered, file = "03_extractedData/190321_KO-RNAseq/DEseq/onlyDifferentiallyEspressedGenes_allTargets.txt", append = TRUE, quote = FALSE, row.names = FALSE, col.names = TRUE)
    
  }
  
}