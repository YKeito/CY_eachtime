#colnames(CY15_Table_1hTFDEGssort)
#head(CY15_Table_1hTFDEGssort)
CY15_Table_intersection <- CY15_Table_20180801[c(CY15_Table_20180801$CY15_1h_normalize_degree != 0 &
                                                   CY15_Table_20180801$CY15_3h_normalize_degree != 0 &
                                                   CY15_Table_20180801$CY15_12h_normalize_degree != 0 &
                                                   CY15_Table_20180801$CY15_24h_normalize_degree != 0),
                                               c("AGI", "TF_family.", colnames(CY15_Table_20180801)[grep("normalize", colnames(CY15_Table_20180801))])]

test <- list(CY15_Table_1hTFDEGssort, CY15_Table_3hTFDEGssort, CY15_Table_12hTFDEGssort, CY15_Table_24hTFDEGssort, CY15_Table_TFDEGsintersection)
title <- c("CY15_1h_top5TF_attibute", "CY15_3h_top5TF_attibute", "CY15_12h_top5TF_attibute", "CY15_24h_top5TF_attibute", "CY15_intersection_top5TF_attibute")
alldata <- c()
i <- 1
for(i in i:length(test)){
  temp <- c("AGI", colnames(test[[i]])[grep("normalize_degree", colnames(test[[i]]))])
  CY15_top5TF <- test[[i]][1:5, temp]
  T_data <- CY15_Table_20180801[match(CY15_top5TF[, "AGI"], CY15_Table_20180801[, "AGI"]), grep("MCLNum", colnames(CY15_Table_20180801))]
  symbol <- allRNASeq$gene_symbol[match(CY15_top5TF[, "AGI"], rownames(allRNASeq))]
  sample <- rep(title[i], times = 5)
  data <- data.frame(CY15_top5TF, T_data, symbol, sample, stringsAsFactors = F)
  alldata <- rbind(alldata, data)
  i <- i+1
}

alldatadupli <- alldata[which(duplicated(alldata$AGI)), ]
T_sampledupli <- alldatadupli$sample
test <- alldata[which(!duplicated(alldata$AGI)), ]
T_sample <- test$sample[match(alldatadupli$AGI, test$AGI)]
T_sample <- paste0(T_sample, T_sampledupli)
test$sample[match(alldatadupli$AGI, test$AGI)] <- T_sample
test <- na.omit(test)

output <- paste0("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/tGRN_attribute/CY15/", "CY15_tGRN_attribute", ".txt")
write.table(test, file = output, append=F, quote = F, sep = "\t", row.names = F)
