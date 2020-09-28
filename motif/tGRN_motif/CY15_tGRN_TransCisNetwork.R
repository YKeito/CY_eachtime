"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/motif/tGRN_motif/CY15_tGRN_TransCisNetwork.R"
Sample <- c(rep("CY15_1h", each = 5*length(unique(names(unlist(allMotif_ID[[1]])[!is.na(unlist(allMotif_ID[[1]]))])))),
            rep("CY15_3h", each = 5*length(unique(names(unlist(allMotif_ID[[2]])[!is.na(unlist(allMotif_ID[[2]]))])))),
            rep("CY15_12h", each = 5*length(unique(names(unlist(allMotif_ID[[3]])[!is.na(unlist(allMotif_ID[[3]]))])))),
            rep("CY15_24h", each = 5*length(unique(names(unlist(allMotif_ID[[4]])[!is.na(unlist(allMotif_ID[[4]]))]))))
)
MCLNum <- c(rep(unique(names(unlist(allMotif_ID[[1]])[!is.na(unlist(allMotif_ID[[1]]))])), each = 5),
            rep(unique(names(unlist(allMotif_ID[[2]])[!is.na(unlist(allMotif_ID[[2]]))])), each = 5),
            rep(unique(names(unlist(allMotif_ID[[3]])[!is.na(unlist(allMotif_ID[[3]]))])), each = 5),
            rep(unique(names(unlist(allMotif_ID[[4]])[!is.na(unlist(allMotif_ID[[4]]))])), each = 5)
)
CY15_summary <- data.frame(CY15_MotifClustering,
                           CY = Sample,
                           MCLNum = MCLNum,
                           stringsAsFactors = F
)
CY15_summary <- CY15_summary[CY15_summary$score > 0, ]

####trans####
condition <- c("CY15_1h", "CY15_3h", "CY15_12h", "CY15_24h")
test <- CY15_Table_20180801[CY15_Table_20180801$TF_family. != "No", c(1, 2, grep("MCL", colnames(CY15_Table_20180801)), grep("normalize", colnames(CY15_Table_20180801)))]
allhub <- list()
n <- 1
for(n in n:length(condition)){
  test <- test[order(test[, paste0(condition[n], "_normalize_degree")], decreasing = T), ]
  total <- unique(test[, paste0(condition[n], "_MCLNum")][test[, paste0(condition[n], "_MCLNum")] != 0][!is.na(test[, paste0(condition[n], "_MCLNum")][test[, paste0(condition[n], "_MCLNum")] != 0])])
  total <- total[order(total)]
  T_hub <- c()
  i <- 1
  for(i in i:length(total)){
    T_hub <- c(T_hub, test$TF_family.[which(total[i] == test[, paste0(condition[n], "_MCLNum")])][1])
    names(T_hub)[i] <- paste0("MCLNum", total[i])
    i <- i+1
  }
  allhub <- c(allhub, list(T_hub))
  n <- n+1
}
names(allhub) <- condition

####Cis-elemnt####
temp <- unique(CY15_summary$Cis)[unique(CY15_summary$Cis) != "Other"]
allCis <- list(CY15_summary$Sample,
               CY15_summary$Sample[CY15_summary$CY != "CY15_1h"],
               CY15_summary$Sample[CY15_summary$CY == "CY15_12h" | CY15_summary$CY == "CY15_24h"],
               CY15_summary$Sample[CY15_summary$CY == "CY15_24h"]
)
names(allCis[[1]]) <- CY15_summary$Cis
names(allCis[[2]]) <- CY15_summary$Cis[CY15_summary$CY != "CY15_1h"]
names(allCis[[3]]) <- CY15_summary$Cis[CY15_summary$CY == "CY15_12h" | CY15_summary$CY == "CY15_24h"]
names(allCis[[4]]) <- CY15_summary$Cis[CY15_summary$CY == "CY15_24h"]

temp <- c("bHLH", "WRKY", "MYB")
time <- c("1h", "3h", "12h", "24h")
network <- c()
interaction_value <- c()
n <- 1
for(n in n:length(temp)){
  source <- list()
  target <- list()
  ttime <- c()
  m <- 1
  for(m in m:length(allhub)){
    Trans <- allhub[[m]][which(allhub[[m]] == temp[n])]
    Trans <- paste0("CY15_", time[m], "_", names(Trans), ",", Trans)
    T_Cis <- allCis[[m]][which(names(allCis[[m]]) == temp[n])]
    T_Cis <- paste0(T_Cis, ",", names(T_Cis))
    Cis <- unique(T_Cis)
    source <- c(source, list(rep(Trans, each = length(Cis))))
    target <- c(target, list(rep(Cis, times = length(Trans))))
    ttime <- c(ttime, paste0("CY15_", rep(time[m], each = length(source[[m]]))))
    m <- m+1
  }
  
  source <- unlist(source)
  target <- unlist(target)
  T_source <- c()
  T_target <- c()
  attr <- c()
  o <- 1
  for(o in o:length(source)){
    itiji <- str_split(source, pattern = ",")[[o]]
    T_source <- c(T_source, itiji[1])
    T_target <- c(T_target, str_split(target, pattern = ",")[[o]][1])
    attr <- c(attr, itiji[2])
  }
  
  network <- rbind(network,   data.frame(Trans = T_source, 
                                         Cis =T_target,
                                         time = ttime,
                                         attr_Trans = attr,
                                         stringsAsFactors = F)
  )
}

network <- network[grep("MCLNum", network$Trans, invert = F), ]
network <- network[grep("MCLNum", network$Cis, invert = F), ]

Node <- union(network$Trans, network$Cis)
NodeColour <- rep("", times = length(Node))
names(NodeColour) <- Node

TransCis <- intersect(network$Trans, network$Cis)
Cis <- setdiff(intersect(Node, network$Cis), TransCis)
Trans <- setdiff(intersect(Node, network$Trans), TransCis)

NodeColour[match(TransCis, names(NodeColour))] <- "TransCis"
NodeColour[match(Cis, names(NodeColour))] <- "Cis"
NodeColour[match(Trans, names(NodeColour))] <- "Trans"

Nodetime <- c()
i <- 1
for(i in i:length(Node)){
  Nodetime <- c(Nodetime, str_split(Node, pattern = "_MCL")[[i]][1])
}
Nodetime[Nodetime == "CY15_1h"] <- "CY15_01h"
Nodetime[Nodetime == "CY15_3h"] <- "CY15_03h"

Node_TF <- rep("", times = length(Node))


CY15_tGRN_attribute <- data.frame(Node = Node,
                                  Nodetime = Nodetime,
                                  NodeColour = NodeColour,
                                  stringsAsFactors = F
)
rownames(CY15_tGRN_attribute) <- c() 

write.table(network, "~/Nakano_RNAseq/network_analysis/result_Table/tGRN_Cis_network/CY15_tGRN_Cis_network.txt", sep = "\t", quote = F, row.names = F)
write.table(CY15_tGRN_attribute, "~/Nakano_RNAseq/network_analysis/result_Table/tGRN_Cis_network/CY15_tGRN_attribute.txt", sep = "\t", quote = F, row.names = F)
