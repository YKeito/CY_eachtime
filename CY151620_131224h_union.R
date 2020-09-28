#"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/CY151620_131224h_union.R"
allRNASeq <- read.table("~/Nakano_RNAseq/network_analysis/base/CY_base/allRNASeq_union.txt", sep = "\t", row.names = 1, header = T, stringsAsFactors = F)
####calculate PCC####
allRNASeq_foldchange <- allRNASeq[, 1:10]

n <- 1
m <- 2
base <- c()
PCC <- c()
PCC_all <- list()
PCC_pvalue <- c()
PCC_pvalue_all <- list()
total <- nrow(allRNASeq_foldchange)
for(n in n:c(total-1)){
  for(m in m:total){
    base <- cor.test(as.numeric(allRNASeq_foldchange[n, ]), as.numeric(allRNASeq_foldchange[m, ]), method = "pearson")
    PCC <- c(PCC, base$estimate)
    PCC_pvalue <- c(PCC_pvalue, base$p.value)
    m <- m+1
  }
  PCC_all <- c(PCC_all, list(PCC))
  PCC <- c()
  PCC_pvalue_all <- c(PCC_pvalue_all, list(PCC_pvalue))
  PCC_pvalue <- c()
  print(total-n)
  n <- n+1
  m <- n+1
}

####PCC q-value####
PCC_qvalue_all <- p.adjust(unlist(PCC_pvalue_all), method = "BH")
#####Cytoscape_format####
source_genes <- combn(rownames(allRNASeq_foldchange), 2)[1, ]
target_genes <- combn(rownames(allRNASeq_foldchange), 2)[2, ]

allRNASeq_foldchange_cytoscape <- data.frame(source_genes = source_genes, 
                                             interaction_value = unlist(PCC_all), 
                                             target_genes = target_genes,
                                             p_value = unlist(PCC_pvalue_all),
                                             q_value = PCC_qvalue_all
)

allRNASeq_cytoscape_th <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.05, ]
allRNASeq_cytoscape_th_possitive <- allRNASeq_cytoscape_th[allRNASeq_cytoscape_th$interaction_value > 0, ]
write.table(allRNASeq_cytoscape_th, file = "~/Nakano_RNAseq/network_analysis/cytoscape/CY151620_131224h_union.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(allRNASeq_cytoscape_th_possitive, file = "bigdata/yasue/tGRN_Groping/sGRN_possitive.txt", append=F, quote = F, sep = "\t", row.names = F)