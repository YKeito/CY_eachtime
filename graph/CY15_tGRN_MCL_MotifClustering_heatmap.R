"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/graph/CY15_tGRN_MCL_MotifClustering_heatmap.R"

Cluster_MotifID <- list()
m <- 1
for(m in m:length(unique(CY15_data$res.hc.cluster))){
  Cluster_MotifID <- c(Cluster_MotifID, list(rownames(CY15_data)[CY15_data$res.hc.cluster == m]))
  m <- m+1
}

i <- 1
total_i <- length(unique(CY15_MotifID_consensus$Name))
allcount <- list()
for(i in i:total_i){
  target <- CY15_MotifID_consensus$MotifID[CY15_MotifID_consensus$Name == unique(CY15_MotifID_consensus$Name)[i]]
  count <- c()
  n <- 1
  for(n in n:length(Cluster_MotifID)){
    count <- c(count, length(intersect(Cluster_MotifID[[n]], target)))
    n <- n+1
  }
  names(count) <- paste0("MotifCluster", unique(CY15_data$res.hc.cluster))
  allcount <- c(allcount, list(count/sum(count)))
  i <- i+1
}
names(allcount) <- unique(CY15_MotifID_consensus$Name)

CY15_MotifClustering <- c()
CY15_MotifClustering <- data.frame(Sample = rep(unique(CY15_MotifID_consensus$Name), each = 5), 
                     score = unlist(allcount),
                     MotifCluster = rep(paste0("MotifCluster", unique(CY15_data$res.hc.cluster)), times = total_i),
                     sample_sort = rep(length(unique(CY15_MotifID_consensus$Name)):1, each = 5),
                     Cis = rep(c("bHLH", "MYB", "Other", "Other", "WRKY"), times = total_i),
                     stringsAsFactors = F
                     )
rownames(CY15_MotifClustering) <- c()

library(ggplot2)
g <- ggplot(CY15_MotifClustering, aes(x = MotifCluster, y = reorder(Sample, sample_sort), fill = score))
g <- g + geom_tile(color = "black", size = 0.1)
g <- g + scale_fill_gradient2(high = "red")
g <- g+theme_dark()
g <- g + theme_linedraw()
#g <- g + coord_flip()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(g)
ggsave("~/Nakano_RNAseq/network_analysis/results_Fig/tGRN_Subcluster_MotifID_Enrichment/CY15/CY15_tGRN_Subcluster_MotifID_Cluster_heatmap.png", g, width = 8, height = 5)