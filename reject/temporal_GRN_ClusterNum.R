#1, primary GRNのTF top5 list
data <- na.omit(CY15_Table_20180801)
primary_GRNTF <- CY15_Table_1hTFDEGssort$AGI[1:5]
primary_GRNMCLNum <- data$CY15_1h_MCLNum[match(primary_GRNTF, data$AGI)]
primary_GRNMCLNumuni <- unique(primary_GRNMCLNum)

allok <- c()
allTF <- list()

ClusterNum1 <- c()
ClusterNum2 <- c()
ClusterNum3 <- c()
ClusterNum4 <- c()
i <- 1
for(i in i:length(primary_GRNMCLNum)){
  temp <- data[data$CY15_1h_MCLNum == primary_GRNMCLNum[i], ] #primart_GRNTFが含まれるClusterのノードを抽出
  temp <- temp[temp$CY15_3h_DEGs != "No" & temp$TF_family. != "No", ] #そのクラスターにおいて3hでTFかつDEGsのnodeをピックアップ
  temp <- temp[order(temp$CY15_3h_normalize_degree, decreasing = T), ]
  if(nrow(temp) != 0){
    ClusterNum1 <- c(ClusterNum1, primary_GRNMCLNum[i])
    secondary_GRNTF <- temp$AGI[1:5] #secondary_GRNのtop5のTFを抽出
    secondary_GRNMCLNum <- data$CY15_3h_MCLNum[match(secondary_GRNTF, data$AGI)] #top5のTFがあるクラスター番号をピックアップ
    }else{
      break
    }
  n <- 1
  for(n in n:length(secondary_GRNMCLNum)){
    temp2 <- data[data$CY15_3h_MCLNum == secondary_GRNMCLNum[n], ] #primart_GRNTFが含まれるClusterのノードを抽出
    temp2 <- temp2[temp2$CY15_12h_DEGs != "No" & temp2$TF_family. != "No", ] #そのクラスターにおいて3hでTFかつDEGsのnodeをピックアップ
    temp2 <- temp2[order(temp2$CY15_12h_normalize_degree, decreasing = T), ]
    if(nrow(temp2) != 0){
      ClusterNum2 <- c(ClusterNum2, secondary_GRNMCLNum[n])
      third_GRNTF <- temp2$AGI[1:5] #fouth_GRNTFの候補TFをピックアップ
      third_GRNMCLNum <- data$CY15_12h_MCLNum[match(third_GRNTF, data$AGI)] #top5が入っているGRNを確認
      }else{
        break
        }
    m <- 1
    for(m in m:length(third_GRNMCLNum)){
      temp3 <- data[data$CY15_12h_MCLNum == third_GRNMCLNum[m], ] #primart_GRNTFが含まれるClusterのノードを抽出
      temp3 <- temp3[temp3$CY15_24h_DEGs != "No" & temp3$TF_family. != "No", ] #そのクラスターにおいて3hでTFかつDEGsのnodeをピックアップ
      temp3 <- temp3[order(temp3$CY15_24h_normalize_degree, decreasing = T), ]
      if(nrow(temp3) != 0){
        ClusterNum3 <- c(ClusterNum3, third_GRNMCLNum[m])
        fouth_GRNTF <- temp3$AGI[1:5]
        fouth_GRNMCLNum <- data$CY15_24h_MCLNum[match(fouth_GRNTF, data$AGI)]
        }else{
          break
          }
      o <- 1
      for(o in o:length(fouth_GRNMCLNum)){
        ClusterNum4 <- c(ClusterNum4, fouth_GRNMCLNum[o])
        o <- o+1
        }
      m <- m+1
    }
    n <- n+1
  }
  print(i)
  i <- i+1
}
