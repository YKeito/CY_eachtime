#"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/graph/tGRN_MCL_AverageOfExpression.R"
CYlist <- list(CY15_Table_20180801, CY16_Table_20180801, CY20_Table_20180801)
CY <- c("CY15", "CY16", "CY20")
n <- 1
for(n in n:length(CYlist)){
  CYData <- paste0(CY[n], c("_1h", "_3h", "_12h", "_24h"))
  allMCLNode <- list()
  allMCLNum <- list()
  AverageOfExpression <- list()
  TF <- list()
  NumTF <- list()
  pvalue <- list()
  
  k <- 1
  for(k in k:length(CYData)){
    temp <- c(CYlist[[n]][, "TF_family."] != "No" & CYlist[[n]][, paste0(CYData[k], "_DEGs")] != "No")
    temp2 <- c(grep("AGI", colnames(CYlist[[n]])), grep(CYData[k], colnames(CYlist[[n]])))
    
    T_data <- CYlist[[n]][temp, temp2]
    T_data <- na.omit(T_data)
    T_data <- T_data[T_data[, paste0(CYData[k], "_MCLNum")] != 0, ]
    
    T_MCLNum <- sort(unique(T_data[, paste0(CYData[k], "_MCLNum")]))
    
    T_MCLNode <- c()
    T_data <- c()
    T_TF <- c()
    T_NumTF <- c()
    T_pvalue <- c()
    i <- 1
    for(i in i:length(T_MCLNum)){
      T_MCLNode <- c(T_MCLNode, sum(CYlist[[n]][, paste0(CYData[k], "_MCLNum")] == T_MCLNum[i], na.rm = T))
      names(T_MCLNode[i]) <- paste0("sub", T_MCLNum[i])
      temp <- CYlist[[n]][, paste0(CYData[k], "_MCLNum")] == T_MCLNum[i]
      temp <- CYlist[[n]][which(temp == "TRUE"), ]
      T_data <- c(T_data, mean(allRNASeq[match(temp$AGI, rownames(allRNASeq)), CYData[k]]))
      names(T_data)[i] <- paste0(CYData[k], "_MCL0", T_MCLNum[i])
      T_NumTF <- c(T_NumTF, sum(temp[, "TF_family."] != "No"))
      T_TF <- c(T_TF, paste0(temp$AGI[temp[, "TF_family."] != "No"], collapse = " | "))
      a <- length(intersect(temp$AGI, temp$AGI[temp[, "TF_family."] != "No"]))
      b <- sum(CYlist[[n]][, "TF_family."] != "No")-a
      c <- sum(temp[, "TF_family."] == "No")
      d <- sum(CYlist[[n]][, "TF_family."] != "No" & CYlist[[n]][, paste0(CYData[k], "_DEGs")] != "No")
      mx <- matrix(c(a, b, c, d), nrow=2, byrow=T)
      T_pvalue <- c(T_pvalue, fisher.test(mx)$p.value)
      i  <- i+1
    }
    allMCLNode <- c(allMCLNode, list(T_MCLNode))
    names(allMCLNode[k]) <- CYData[k]
    AverageOfExpression <- c(AverageOfExpression, list(T_data))
    names(AverageOfExpression[k]) <- CYData[k]
    TF <- c(TF, list(T_TF))
    names(TF[k]) <- CYData[k]
    NumTF <- c(NumTF, list(T_NumTF))
    names(NumTF[k]) <- CYData[k]
    pvalue <- c(pvalue, list(T_pvalue))
    
    m <- 1
    T_all <- c()
    for(m in m:length(T_MCLNum)){
      temp <- "0000"
      temp <- paste0(substr(temp, 1, nchar(temp)-nchar(T_MCLNum[m])), T_MCLNum[m])
      T_all <- c(T_all, temp)
    }
    allMCLNum <- c(allMCLNum, list(T_all))
    names(allMCLNum[k]) <- CYData[k]
    print(k)
    k <- k+1
  }
  
  CYData <- c("CY15_01h", "CY15_03h", "CY15_12h", "CY15_24h")
  data <- data.frame(MCLNum = unlist(allMCLNum),
                     MCLClu = as.character(unlist(allMCLNode)),
                     sample = c(rep(CYData[1], times = length(unlist(allMCLNum[[1]]))),
                                rep(CYData[2], times = length(unlist(allMCLNum[[2]]))),
                                rep(CYData[3], times = length(unlist(allMCLNum[[3]]))),
                                rep(CYData[4], times = length(unlist(allMCLNum[[4]])))
                     ),
                     AverageOfExpression = unlist(AverageOfExpression),
                     TF = unlist(TF),
                     NumTF = unlist(NumTF),
                     pvalue = unlist(pvalue),
                     log2_pvalue = -log2(unlist(pvalue)),
                     stringsAsFactors = F
  )
  
  temp <- unique(data$sample)
  m <- 1
  for(m in m:length(temp)){
    T_data <- data[data$sample == temp[m], ]
    library(ggplot2)
    #heatmap
    g <- ggplot(T_data, aes(x = MCLNum, y = sample, fill = AverageOfExpression))
    g <- g + geom_tile()
    g <- g + theme_bw()
    g <- g + scale_fill_gradient2(low = "blue", high = "red")
    g <- g + coord_flip()
    plot(g)
    title <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/tGRN_MCL_AverageEx/", CY[n], "/", temp[m], "_tGRN_AverageEx.png")
    ggsave(file = title, plot = g)
  }
  title <- paste0("~/Nakano_RNAseq/network_analysis/eachCY_Table/temporal_network_tGRN/", CY[n], "/", CY[n], "_tGRN.txt")
  write.table(data, file = title, append=F, quote = F, sep = "\t", row.names = F)
  assign(paste0(CY[n], "_tGRN_MCL"), data)
  n <- n+1
}