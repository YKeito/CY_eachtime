#"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/motif/CYMCLNumDEGs_ame_sequences_pvalue.R"
#library("stringr")
all <- unlist(allMotif_ID)
all <- all[!is.na(all)]
time <- c("1h", "3h", "12h", "24h")
T_pvalue <- c()
obname <- c()
Motif_name <- c()
Cluster <- c()
T_time <- c()
m <- 1
for(m in m:length(allMotif_ID)){
  temp <- allMotif_ID[[m]][!is.na(allMotif_ID[[m]])]
  target <- unique(names(temp))
  n <- 1
  for(n in n:length(target)){
    i <- 1
    for(i in i:length(unique(all))){
      nn <- length(temp[names(temp) == target[n]])
      a <- length(grep(unique(all)[i], temp[names(temp) == target[n]], invert = F))
      MM <- sum(all == unique(all)[i], na.rm = T)
      b <- MM-a
      c <- nn-a
      d <- length(all)-a-b-c
      NN <-length(all)
      mx <- matrix(c(a, b, c, d), nrow=2, byrow=T)
      T_pvalue <- c(T_pvalue, fisher.test(mx, alternative = "greater")$p.value)
      obname <- c(obname, c("CY15"))
      Cluster <- c(Cluster, target[n])
      T_time <- c(T_time, time[m])
      Motif_name <- c(Motif_name, unique(all)[i])
      i <- i+1
    }
    n <- n+1
  }
  print(m)
  m <- m+1
}

temp <- paste0(obname, "_", T_time, "_", Cluster)
sort <- c()
m <- 1
for(m in m:length(unique(temp))){
  sort <- c(sort, rep(m, times = sum(unique(temp)[m] == temp)))
}

data <- data.frame(sort = sort,
                   CY_MCL = paste0(obname, "_", T_time, "_", Cluster),
                   Motif_ID = Motif_name,
                   pvalue = T_pvalue,
                   enrichment = -log10(T_pvalue),
                   stringsAsFactors = F
                   )

data <- data[data$pvalue < 0.05, ]

g <- ggplot(data, aes(x = Motif_ID, y = reorder(CY_MCL, sort), fill = enrichment))
g <- g + geom_tile()
#g <- g + theme_bw()
g <- g + scale_fill_gradient2(high = "red")
#g <- g + coord_flip()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(g)


ggsave("~/Nakano_RNAseq/network_analysis/results_Fig/tGRN_Subcluster_MotifID_Enrichment/CY15_tGRN_Subcluster_MotifID_Enrichment_heatmap.png", g)
write.table(data, file = "~/Nakano_RNAseq/network_analysis/result_Table/tGRN_Subcluster_MotifID/tGRN_Subcluster_MotifID_summary.txt", append=F, quote = F, sep = "\t", row.names = F)