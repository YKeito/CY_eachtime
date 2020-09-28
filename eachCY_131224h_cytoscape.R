mastercluster_pair <- paste0(allRNASeq_cytoscape_th$source_genes, allRNASeq_cytoscape_th$target_genes)
mastercluster_pair <- data.frame(AGI_pair = mastercluster_pair,
                                 PCC = allRNASeq_cytoscape_th$interaction_value
                                 )

AGI_list <- list(rownames(CY15_1h_FDR005), rownames(CY15_3h_FDR005), rownames(CY15_12h_FDR005), rownames(CY15_24h_FDR005),
                  rownames(CY16_1h_FDR005), rownames(CY16_3h_FDR005), rownames(CY16_12h_FDR005), rownames(CY16_24h_FDR005),
                  rownames(CY20_1h_FDR005), rownames(CY20_3h_FDR005), rownames(CY20_12h_FDR005), rownames(CY20_24h_FDR005)
                  )
names(AGI) <- list("CY15_1h_FDR005", "CY15_3h_FDR005", "CY15_12h_FDR005", "CY15_24h_FDR005",
                   "CY16_1h_FDR005", "CY16_3h_FDR005", "CY16_12h_FDR005", "CY16_24h_FDR005",
                   "CY20_1h_FDR005", "CY20_3h_FDR005", "CY20_12h_FDR005", "CY20_24h_FDR005"
                   )

i <- 1
for(i in i:length(AGI_list)){
  CY <- combn(AGI_list[[i]], 2)
  pair <- paste0(CY[1, ], CY[2, ])
  edge <- intersect(mastercluster_pair$AGI_pair, pair)
  PCC <- mastercluster_pair$PCC[match(edge, mastercluster_pair$AGI_pair)]
  souce <- c()
  target <- c()
  m <- 1
  for(m in m:length(edge)){
    souce <- c(souce, substr(edge[m], 1, 9))
    target <- c(target, substr(edge[m], 10, 18))
    m <- m+1
  }
  temp <- data.frame(souce = souce,
                     PCC = PCC,
                     target = target
                     )
  File <- paste0("~/Nakano_RNAseq/network_analysis/base/eachCY_131224h/", names(AGI)[i])
  write.table(temp, file = File, append=F, quote = F, sep = "\t", row.names = F)
  print(i)
}
