#"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/graph/smoothgraph/CY_temporal_smoothgraph.R"
library(ggplot2)
library(reshape2)
Num <- c(7)
condition <- c("CY15_1h_MCLNum")
m <- 1
for(m in m:length(Num)){
  n <- 1
  for(n in n:length(condition)){
    T_AGI <- na.omit(CY15_Table_20180801$AGI[CY15_Table_20180801[, condition[n]] == Num[m]])
    temp <- str_split(condition[n], pattern = "_MCL")[[1]][1]
    T_data <- allRNASeq[match(T_AGI, rownames(allRNASeq)), ]
    T_data <- T_data[order(T_data[, temp], decreasing = T), ]
    sort <- 1:nrow(T_data)
    T_expression <- c(rep(0, times = length(T_data$CY15_1h)),
                      T_data$CY15_1h,
                      T_data$CY15_3h,
                      T_data$CY15_12h,
                      T_data$CY15_24h
    )
    df <- data.frame(AGI = rep(rownames(T_data), time = 5),
                     expression = T_expression,
                     Category = rep(c("0h", "01h", "03h", "12h", "24h"), each = length(T_AGI)),
                     time = rep(c(0, 1, 3, 12, 24), each = length(T_AGI)),
                     group = rep("", times = length(T_expression)),
                     AGI_sort = rep(sort, time = 5),
                     stringsAsFactors = F
    )
    i <- 1
    PCC <- c()
    for(i in i:length(unique(df$AGI))){
      PCC <- c(PCC, cor(df$expression[df$AGI == unique(df$AGI)[1]], df$expression[df$AGI == unique(df$AGI)[i]]))
    }
    names(PCC) <- unique(df$AGI)
    i <- 1
    for(i in i:length(PCC)){
      df$group[grep(names(PCC)[PCC < 0][i], df$AGI)] <- paste0("02negative", "_", sum(PCC < 0), " genes")
      df$group[grep(names(PCC)[PCC > 0][i], df$AGI)] <- paste0("01possitive", "_", sum(PCC > 0), " genes")
    }
    df$expression[df$expression < -2] <- -2
    df$expression[df$expression > 2] <- 2
    g1 <- ggplot(df, aes(x = Category, y = reorder(AGI, AGI_sort), fill = expression))
    g1 <- g1 + geom_tile()
    g1 <- g1 + theme_bw()
    g1 <- g1 + scale_fill_gradient2(limits = c(-2,2), low = "blue", high = "red", na.value = "white")
    g1 <- g1 + ggtitle(paste0(condition[n], Num[m]))
    plot(g1)
    title <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/tGRN_analysis/expression_graph/heatmap/", condition[n], Num[m], "_heatmap", ".png")
    ggsave(title, g1)
    
    ####geom_smooth####
    g <- ggplot(data = df, aes(x = time, y = expression, colour = group))
    #g <- g + geom_point(mapping = aes(x = time, y = expression, colour = group))
    #g <- g + geom_line(mapping = aes(x = time, y = expression, group = AGI))
    g <- g + geom_smooth(method = "loess", mapping = aes(x = time, y = expression, colour = group), level = 0.95)
    g <- g + theme(legend.position="top")
    g <- g + ggtitle(paste0(condition[n], Num[m]))
    plot(g)
    title <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/tGRN_analysis/expression_graph/smoothgraph/", condition[n], Num[m], "_smooth", ".png")
    ggsave(title, g)
    
    AGI_correlation <- dcast(df, AGI ~ group, value.var = "group")
    title <- paste0("~/Nakano_RNAseq/network_analysis/result_Table/tGRN_MCL_expression/", condition[n], Num[m], "AGI_correlation.txt")
    write.table(AGI_correlation, title, sep = "\t", quote = F, row.names = F)
  }
}
