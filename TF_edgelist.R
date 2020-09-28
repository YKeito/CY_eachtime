#CY151620_131224h_cytoscape_th <- read.table(file = "~/Nakano_RNAseq/network_analysis/base/eachCY_131224h/CY151620_131224h_unionPCC.txt", sep = "\t", stringsAsFactors = F, header = T)
#AGI <- union(CY151620_131224h_cytoscape_th$source_genes, CY151620_131224h_cytoscape_th$target_genes)
#check
#setdiff(TF_AGI, AGI)
#setdiff(rownames(allRNASeq_foldchange), union(CY151620_131224h_cytoscape_th$source_genes, CY151620_131224h_cytoscape_th$target_genes))
#

#object <- list(CY15_Table_1hTFDEGssort, CY15_Table_3hTFDEGssort, CY15_Table_12hTFDEGssort, CY15_Table_24hTFDEGssort, CY15_Table_TFDEGsintersection)

TF_all <- rep(0, times = nrow(CY151620_131224h_cytoscape_th))
temp <- intersect(TF_family$AGI, CY151620_131224h_cytoscape_th$source_genes)
TF_all[match(temp, CY151620_131224h_cytoscape_th$source_genes)] <- TF_family$TF[match(temp, TF_family$AGI)]

temp <- intersect(TF_family$AGI, CY151620_131224h_cytoscape_th$target_genes)
TF_all[match(temp, CY151620_131224h_cytoscape_th$target_genes)] <- TF_family$TF[match(temp, TF_family$AGI)]

File <- c("CY15_Table_1hTFDEGs_1hsort", "CY15_Table_3hTFDEGs_3hsort", "CY15_Table_12hTFDEGs_12hsort", "CY15_Table_24hTFDEGs_24hsort", "CY15_Table_TFDEGsintersection_1hsort")
test <- list(CY15_Table_1hTFDEGssort$AGI, CY15_Table_3hTFDEGssort$AGI, CY15_Table_12hTFDEGssort$AGI, CY15_Table_24hTFDEGssort$AGI, CY15_Table_TFDEGsintersection$AGI)
i <- 1
for(i in i:length(test)){
  TF_AGI <- test[[i]]
  
  souce_genes_TF <- rep(0, times = nrow(CY151620_131224h_cytoscape_th))
  target_genes_TF <- rep(0, times = nrow(CY151620_131224h_cytoscape_th))
  
  total <- length(TF_AGI)
  n <- 1
  for(n in n:total){
    temp <- grep(TF_AGI[n], CY151620_131224h_cytoscape_th$source_genes)
    souce_genes_TF[temp] <- "TF"
    
    temp <- grep(TF_AGI[n], CY151620_131224h_cytoscape_th$target_genes)
    target_genes_TF[temp] <- "TF"
    print(n)
    n <- n+1
  }
  
  TF_target <- data.frame(CY151620_131224h_cytoscape_th, 
                          souce_genes_TF,
                          target_genes_TF,
                          TF_all = TF_all,
                          stringsAsFactors = F
  )
  TF_target <- TF_target[, c(1, 6, 2, 4, 5, 3, 7, 8 )]
  
  TF_target <- TF_target[TF_target$souce_genes_TF != 0 | TF_target$target_genes_TF != 0, ]
  output <- paste0(File[i], ".txt")
  
  write.table(TF_target, file=paste0("~/Nakano_RNAseq/network_analysis/eachCY_Table/CY15_TF_DGEs_sort/TF_edgelist/", output), append=F, quote=F, sep="\t", row.names = F, col.names = T)
  print(i)
  i <- i+1
}
