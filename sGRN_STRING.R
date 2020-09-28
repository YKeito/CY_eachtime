#STRING <- read.table("~/bigdata/yasue/STRING/3702.protein.links.full.v10.5.txt", stringsAsFactors = F, header = T)
#check
#sum(nchar(STRING$protein1) != 16)
#sum(nchar(STRING$protein2) != 16)
#STRING$protein1 <- substr(STRING$protein1, 6, 14)
#STRING$protein2 <- substr(STRING$protein2, 6, 14)
#STRING <- STRING[, c("protein1", "protein2", "combined_score")]


CYpair <- paste0(CY151620_131224h_cytoscape_th$source_genes, CY151620_131224h_cytoscape_th$target_genes)
CY_STRING <- rep("No", times = nrow(CY151620_131224h_cytoscape_th))
CY_STRING[match(intersect(STRINGpair, CYpair), CYpair)] <- STRING$combined_score[match(intersect(STRINGpair, CYpair), STRINGpair)]

TF_source <- rep("No", times = nrow(CY151620_131224h_cytoscape_th))
TF_target <- rep("No", times = nrow(CY151620_131224h_cytoscape_th))
i <- 1
for(i in i:length(TF_family$AGI)){
  TF_source[grep(TF_family$AGI[i], CY151620_131224h_cytoscape_th$source_genes)] <- "Yes"
  TF_target[grep(TF_family$AGI[i], CY151620_131224h_cytoscape_th$target_genes)] <- "Yes"
  print(i)
  i <- i+1
}
data <- data.frame(CY151620_131224h_cytoscape_th[, c(1:3)],
                   TF_source,
                   TF_target,
                   STRING_combined_score = CY_STRING,
                   stringsAsFactors = F)
write.table(data, file = "~/Nakano_RNAseq/network_analysis/add_STRING/sGRN_CY151620PCC.txt", sep = "\t", , append=F, quote = F, row.names = F)

