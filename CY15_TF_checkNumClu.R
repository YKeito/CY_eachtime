CYData <- c("CY15_1h", "CY15_3h", "CY15_12h", "CY15_24h")
allMCLNode <- list()
allMCLNum <- list()
n <- 1
for(n in n:length(CYData)){
  temp <- c(CY15_Table_20180801$TF_family. != "No" & CY15_Table_20180801[, paste0(CYData[n], "_DEGs")] != "No")
  temp2 <- c(grep("AGI", colnames(CY15_Table_20180801)), grep(CYData[n], colnames(CY15_Table_20180801)))
  
  T_data <- CY15_Table_20180801[temp, temp2]
  T_data <- na.omit(T_data)
  T_data <- T_data[T_data[, paste0(CYData[n], "_MCLNum")] != 0, ]
  
  T_MCLNum <- sort(unique(T_data[, paste0(CYData[n], "_MCLNum")]))
  
  T_MCLNode <- c()
  i <- 1
  for(i in i:length(T_MCLNum)){
    T_MCLNode <- c(T_MCLNode, sum(CY15_Table_20180801[, paste0(CYData[n], "_MCLNum")] == T_MCLNum[i], na.rm = T))
    names(T_MCLNode[i]) <- paste0("sub", T_MCLNum[i])
    i  <- i+1
  }
  allMCLNode <- c(allMCLNode, list(T_MCLNode))
  names(allMCLNode[n]) <- CYData[n]
  allMCLNum <- c(allMCLNum, list(T_MCLNum))
  names(allMCLNum[n]) <- CYData[n]
  print(n)
}


data <- data.frame(MCLNum = paste0("MCLNum", unlist(allMCLNum)),
                   MCLClu = unlist(allMCLNode),
                   sample = c(rep(CYData[1], times = length(unlist(allMCLNum[[1]]))),
                              rep(CYData[2], times = length(unlist(allMCLNum[[2]]))),
                              rep(CYData[3], times = length(unlist(allMCLNum[[3]]))),
                              rep(CYData[4], times = length(unlist(allMCLNum[[4]])))
                   )
)
#write.table(data, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/temporal_network_tGRN/CY15/CY15_tGRN.txt", append=F, quote = F, sep = "\t", row.names = F)