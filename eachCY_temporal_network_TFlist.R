#1, primary GRNのTF top5 list
secondary_GRNTF <- c()
secondary_GRNMCLNum <- c()
secondary_normalize_degree <- c()
tertiary_GRNTF <- c()
tertiary_GRNMCLNum <- c()
tertiary_normalize_degree <- c()
quaternary_GRNTF <- c()
quaternary_GRNMCLNum <- c()
quaternary_normalize_degree <- c()

data <- na.omit(CY15_Table_20180801)
T_primary_GRNTF <- CY15_Table_1hTFDEGssort[c(1:5), "AGI"]
T_primary_GRNMCLNum <- data$CY15_1h_MCLNum[match(T_primary_GRNTF, data$AGI)]
T_primary_normalize_degree <- CY15_Table_1hTFDEGssort[c(1:5), "CY15_1h_normalize_degree"]

primary_GRNTF <- rep(T_primary_GRNTF, each = 125)
primary_GRNMCLNum <- rep(T_primary_GRNMCLNum, each =125)
primary_normalize_degree <- rep(T_primary_normalize_degree, each = 125)

n <- 1
for(n in n:length(T_primary_GRNMCLNum)){
  T_secondary_GRNNode <- data[data$CY15_1h_DEGs != "No" & data$TF_family. != "No", ]
  T_secondary_GRN <- T_secondary_GRNNode[T_secondary_GRNNode$CY15_3h_DEGs != "No" & T_secondary_GRNNode$TF_family. != "No", ]
  T_secondary_GRN <- T_secondary_GRN[order(T_secondary_GRN$CY15_3h_normalize_degree, decreasing = T), ]
  T_secondary_GRNTF <- T_secondary_GRN[c(1:5), "AGI"]
  T_secondary_GRNMCLNum <- T_secondary_GRN[c(1:5), "CY15_3h_MCLNum"]
  T_secondary_normalize_degree <- T_secondary_GRN[c(1:5), "CY15_3h_normalize_degree"]
  m <- 1
  for(m in m:length(T_secondary_GRNMCLNum)){
    T_tertiary_GRNNode <- data[data$CY15_3h_DEGs != "No" & data$TF_family. != "No", ]
    T_tertiary_GRN <- T_tertiary_GRNNode[T_tertiary_GRNNode$CY15_12h_DEGs != "No" & T_tertiary_GRNNode$TF_family. != "No", ]
    T_tertiary_GRN <- T_tertiary_GRN[order(T_tertiary_GRN$CY15_12h_normalize_degree, decreasing = T), ]
    T_tertiary_GRNTF <- T_tertiary_GRN[c(1:5), "AGI"]
    T_tertiary_GRNMCLNum <- T_tertiary_GRN[c(1:5), "CY15_12h_MCLNum"]
    T_tertiary_normalize_degree <- T_tertiary_GRN[c(1:5), "CY15_12h_normalize_degree"]
    
    i <- 1
    for(i in i:length(T_tertiary_GRNMCLNum)){
      T_quaternary_GRNNode <- data[data$CY15_12h_DEGs != "No" & data$TF_family. != "No", ]
      T_quaternary_GRN <- T_quaternary_GRNNode[T_quaternary_GRNNode$CY15_24h_DEGs != "No" & T_quaternary_GRNNode$TF_family. != "No", ]
      T_quaternary_GRN <- T_quaternary_GRN[order(T_quaternary_GRN$CY15_24h_normalize_degree, decreasing = T), ]
      T_quaternary_GRNTF <- T_quaternary_GRN[c(1:5), "AGI"]
      T_quaternary_GRNMCLNum <- T_quaternary_GRN[c(1:5), "CY15_24h_MCLNum"]
      T_quaternary_normalize_degree <- T_quaternary_GRN[c(1:5), "CY15_24h_normalize_degree"]
      
      quaternary_GRNTF <- c(quaternary_GRNTF, T_quaternary_GRNTF)
      quaternary_GRNMCLNum <- c(quaternary_GRNMCLNum, T_quaternary_GRNMCLNum)
      quaternary_normalize_degree <- c(quaternary_normalize_degree, T_quaternary_normalize_degree)
      i <- i+1
    }
    tertiary_GRNTF <- c(tertiary_GRNTF, rep(T_tertiary_GRNTF, each = 5))
    tertiary_GRNMCLNum <- c(tertiary_GRNMCLNum, rep(T_tertiary_GRNMCLNum, each = 5))
    tertiary_normalize_degree <- c(tertiary_normalize_degree, T_tertiary_normalize_degree)
    m <- m+1
  }
  secondary_GRNTF <- c(secondary_GRNTF, rep(T_secondary_GRNTF, each = 25))
  secondary_GRNMCLNum <- c(secondary_GRNMCLNum, rep(T_secondary_GRNMCLNum, each =25))
  secondary_normalize_degree <- c(secondary_normalize_degree, T_secondary_normalize_degree)
  print(n)
  n <- n+1
}

T_data <- data.frame(primary_GRNTF,
                     primary_GRNMCLNum,
                     primary_normalize_degree,
                     secondary_GRNTF,
                     secondary_GRNMCLNum,
                     secondary_normalize_degree,
                     tertiary_GRNTF,
                     tertiary_GRNMCLNum,
                     tertiary_normalize_degree,
                     quaternary_GRNTF,
                     quaternary_GRNMCLNum,
                     quaternary_normalize_degree,
                     stringsAsFactors = F
)
T_data <- na.omit(T_data)

write.table(T_data, "~/Nakano_RNAseq/network_analysis/eachCY_Table/temporal_network_TFlist/CY15/CY15_temporalnetwork_TFlist.txt", sep = "\t", , append=F, quote = F, row.names = F)