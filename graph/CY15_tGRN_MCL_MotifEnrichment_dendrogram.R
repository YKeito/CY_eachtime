"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/graph/CY15_tGRN_MCL_MotifEnrichment_dendrogram.R"

library(reshape2)
CY15_Motif_Table <- dcast(df, sample ~ MotifID, value.var = "value")
#http://langstat.hatenablog.com/entry/20141007/1412607600
#install.packages("stringdist")
library(stringdist)
consensus <- CY15_MotifID_consensus$consensu[match(unique(CY15_MotifID_consensus$MotifID), CY15_MotifID_consensus$MotifID)]
names(consensus) <- unique(CY15_MotifID_consensus$MotifID)


LV_Table <- c()
n <- 1
for(n in n:length(consensus)){
  test <- c()
  m <- 1
  for(m in m:length(consensus)){
    test <- c(test, stringdist(consensus[n], consensus[m],  method = "lcs"))
    m <- m+1
  }
  LV_Table <- rbind(LV_Table, test)
  n <- n+1
  m <- m+1
}

rownames(LV_Table) <- names(consensus)
colnames(LV_Table) <- names(consensus)

#http://www.sthda.com/english/articles/28-hierarchical-clustering-essentials/92-visualizing-dendrograms-ultimate-guide/
#http://documentation1458.rssing.com/chan-49683810/all_p2.html
#install.packages("factoextra")
#install.packages("dendextend")
#library(factoextra)
#library(ggplot2)
#library(igraph)

#k:クラスター数の指定
#k_colors, palette:色指定
#show_labels:横軸の名前を表示するか
#color_labels_by_k:クラスターのグループに応じて色を付けるか
#horiz縦軸と横軸を反転するかしないか
#rect:クラスターを囲む
#rect_fill:rectで囲んだ範囲を色塗する
#type:図の形式の指定

res.hc <- eclust(x = LV_Table, 
                 "hclust", 
                 k = 5,
                 method = "euclidean", 
                 graph = FALSE
                 )

g <- fviz_dend(res.hc,
               cex = 1.2,
               color_labels_by_k = TRUE,
               show_labels = TRUE,
               ggtheme = theme_bw(),
               horiz = FALSE,
               rect = TRUE,
               rect_fill = TRUE,
               type = "rectangle"
               )


#fviz_silhouette(res.hc)
# Silhouette width of observation
#sil <- res.hc$silinfo$widths[, 1:3]
# Objects with negative silhouette
#neg_sil_index <- which(sil[, 'sil_width'] < 0)
#sil[neg_sil_index, , drop = FALSE]
#ggsave("~/Nakano_RNAseq/network_analysis/results_Fig/tGRN_Subcluster_MotifID_Enrichment/CY15/CY15_tGRN_MotifID_dendrogram_clustering5.png", g, width = 20, height = 10)

CY15_data <- data.frame(res.hc$cluster)
consensus <- consensus[match(names(consensus), rownames(CY15_data))]
count <- c()
i <- 1
for(i in i:length(unique(CY15_data$res.hc.cluster))){
  count <- c(count, sum(CY15_data$res.hc.cluster == i))
}
print(count)
CY15_data <- data.frame(CY15_data, consensus)
write.table(CY15_data, "~/Nakano_RNAseq/network_analysis/result_Table/tGRN_Subcluster_MotifID/CY15/CY15_tGRN_Subcluster5_MotifID_summary.txt", sep = "\t", quote = F, col.names = T)
