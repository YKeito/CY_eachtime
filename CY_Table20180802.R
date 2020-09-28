#################################################################################################set##########################################################################################
allgenes <- read.table("~/Nakano_RNAseq/network_analysis/base/allgenes_attribute.txt", sep = "\t", header = T, stringsAsFactors = F)
all <- allgenes$AGI
filename <- list.files("~/Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY16", pattern=".txt", full.names = T)
temp2 <- c("12h", "1h", "24h", "3h")
object_name <- c()
object_name_all <- c()
CY16_DEGs <- c()
CY16_DEGs_all <- c()
n <- 1
for (n in 1:length(filename)){
  object_name <- paste0("CY16", "_", temp2[n], "_DEGs")
  object_name_all <- c(object_name_all, object_name)
  object_name <- assign(object_name, read.table(filename[n],  header=T, sep="\t", stringsAsFactors = F))
  CY16_DEGs <- rep("No", times = length(all))
  temp <- intersect(rownames(object_name), all)
  temp <- match(temp, all)
  CY16_DEGs[temp] <- "Yes"
  CY16_DEGs_all <- cbind(CY16_DEGs_all, CY16_DEGs)
  n <- n+1
}
colnames(CY16_DEGs_all) <- object_name_all

#######MCLNum#############################################################################################################################################################################
#cluster info
filename <- list.files("~/Nakano_RNAseq/network_analysis/base/each_CY131224h_masterclusterinfo/CY16", pattern=".csv", full.names = T)
temp2 <- c("12h", "1h", "24h", "3h")
object_name <- c()
object_name_all <- c()
MCL_Num <- c()
MCL_Num_all <- c()
BC <- c()
BC_all <- c()
degree <- c()
degree_all <- c()
normalize_degree <- c()
normalize_degree_all <- c()
n <- 1
for (n in 1:length(filename)){
  object_name <- paste0("CY16", "_", temp2[n], "_MCLNum")
  object_name <- assign(object_name, read.table(filename[n],  header=T, sep=",", stringsAsFactors = F))
  temp1 <- intersect(object_name[, "name"], all)
  temp1 <- match(temp1, all)  
  
  MCL_Num <- rep(0, times = length(all))
  MCL_Num[temp1] <- object_name[, "X__mclCluster"]
  MCL_Num_all <- cbind(MCL_Num_all, MCL_Num)
  
  BC <- rep(0, times = length(all))
  BC[temp1] <- object_name[, "BetweennessCentrality"]
  BC_all <- cbind(BC_all, BC)
  
  degree <- rep(0, times = length(all))
  normalize_degree <- rep(0, times = length(all))
  degree[temp1] <- object_name[, "Degree"]
  normalize_degree[temp1] <- object_name[, "Degree"]/max(object_name[, "Degree"])
  normalize_degree_all <- cbind(normalize_degree_all, normalize_degree)
  degree_all <- cbind(degree_all, degree)
  n <- n+1
}
colnames(MCL_Num_all) <- paste0("CY16_", temp2, "_MCLNum")
colnames(BC_all) <- paste0("CY16_", temp2, "_BetweennessCentrality")
colnames(degree_all) <- paste0("CY16_", temp2, "_degree")
colnames(normalize_degree_all) <- paste0("CY16_", temp2, "_normalize_degree")
##########################################################################################
CY16_Table_20180801 <- data.frame(allgenes, 
                                  CY16_DEGs_all,
                                  MCL_Num_all,
                                  BC_all,
                                  degree_all,
                                  normalize_degree_all,
                                  stringsAsFactors = F
                                  )

CY16_Table_20180801 <- CY16_Table_20180801[, c(1, 6, 2:5, 
                                               grep("1h", colnames(CY16_Table_20180801)), 
                                               grep("3h", colnames(CY16_Table_20180801)),
                                               grep("12h", colnames(CY16_Table_20180801)),
                                               grep("24h", colnames(CY16_Table_20180801))
                                               )]


#write.table(CY16_Table_20180801, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/20180803CY_Table/CY16_Table_20180803.txt", append=F, quote = F, sep = "\t", row.names = F)

#1h
CY16_Table_1hsort <- CY16_Table_20180801[order(CY16_Table_20180801$CY16_1h_normalize_degree, decreasing = T), 
                                         c("AGI", "TF_family.", "CY16_1h_DEGs", colnames(CY16_Table_20180801)[grep("normalize", colnames(CY16_Table_20180801))])]
CY16_Table_1hDEGssort <- CY16_Table_1hsort[CY16_Table_1hsort$CY16_1h_DEGs != "No", ]
CY16_Table_1hTFDEGssort <- CY16_Table_1hDEGssort[CY16_Table_1hDEGssort$TF_family. != "No", ]

#3h
CY16_Table_3hsort <- CY16_Table_20180801[order(CY16_Table_20180801$CY16_3h_normalize_degree, decreasing = T), 
                                         c("AGI", "TF_family.", "CY16_3h_DEGs", colnames(CY16_Table_20180801)[grep("normalize", colnames(CY16_Table_20180801))])]
CY16_Table_3hDEGssort <- CY16_Table_3hsort[CY16_Table_3hsort$CY16_3h_DEGs != "No", ]
CY16_Table_3hTFDEGssort <- CY16_Table_3hDEGssort[CY16_Table_3hDEGssort$TF_family. != "No", ]


#12h
CY16_Table_12hsort <- CY16_Table_20180801[order(CY16_Table_20180801$CY16_12h_normalize_degree, decreasing = T), 
                                         c("AGI", "TF_family.", "CY16_12h_DEGs", colnames(CY16_Table_20180801)[grep("normalize", colnames(CY16_Table_20180801))])]
CY16_Table_12hDEGssort <- CY16_Table_12hsort[CY16_Table_12hsort$CY16_12h_DEGs != "No", ]
CY16_Table_12hTFDEGssort <- CY16_Table_12hDEGssort[CY16_Table_12hDEGssort$TF_family. != "No", ]

#24h
CY16_Table_24hsort <- CY16_Table_20180801[order(CY16_Table_20180801$CY16_24h_normalize_degree, decreasing = T), 
                                         c("AGI", "TF_family.", "CY16_24h_DEGs", colnames(CY16_Table_20180801)[grep("normalize", colnames(CY16_Table_20180801))])]
CY16_Table_24hDEGssort <- CY16_Table_24hsort[CY16_Table_24hsort$CY16_24h_DEGs != "No", ]
CY16_Table_24hTFDEGssort <- CY16_Table_24hDEGssort[CY16_Table_24hDEGssort$TF_family. != "No", ]


CY16_Table_intersection <- CY16_Table_20180801[c(CY16_Table_20180801$CY16_1h_normalize_degree != 0 &
                                                 CY16_Table_20180801$CY16_3h_normalize_degree != 0 &
                                                 CY16_Table_20180801$CY16_12h_normalize_degree != 0 &
                                                 CY16_Table_20180801$CY16_24h_normalize_degree != 0),
                                               c("AGI", "TF_family.", colnames(CY16_Table_20180801)[grep("normalize", colnames(CY16_Table_20180801))])]

CY16_Table_intersection1hsort <- CY16_Table_intersection[order(CY16_Table_intersection$CY16_1h_normalize_degree, decreasing = T), ]
CY16_Table_TFDEGsintersection <- CY16_Table_intersection[CY16_Table_intersection$TF_family. != "No", ]
CY16_Table_TFDEGsintersection <- CY16_Table_TFDEGsintersection[order(CY16_Table_TFDEGsintersection$CY16_1h_normalize_degree, decreasing = T), ]

write.table(CY16_Table_1hTFDEGssort, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/CY16_TF_DGEs_sort/CY16_Table_1hTFDEGs_1hsort.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(CY16_Table_3hTFDEGssort, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/CY16_TF_DGEs_sort/CY16_Table_3hTFDEGs_3hsort.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(CY16_Table_12hTFDEGssort, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/CY16_TF_DGEs_sort/CY16_Table_12hTFDEGs_12hsort.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(CY16_Table_24hTFDEGssort, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/CY16_TF_DGEs_sort/CY16_Table_24hTFDEGs_24hsort.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(CY16_Table_TFDEGsintersection, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/CY16_TF_DGEs_sort/CY16_Table_TFDEGs_intersection_1hsort.txt", append=F, quote = F, sep = "\t", row.names = F)