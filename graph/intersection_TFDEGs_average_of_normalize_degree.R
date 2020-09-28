library(ggplot2)
library(reshape2)
CY_normalizedegree_aver <- data.frame(sample = rep(c("CY15", "CY16", "CY20"), each = 4),
                                      average = c(mean(CY15_Table_TFDEGsintersection$CY15_1h_normalize_degree),
                                                  mean(CY15_Table_TFDEGsintersection$CY15_3h_normalize_degree),
                                                  mean(CY15_Table_TFDEGsintersection$CY15_12h_normalize_degree),
                                                  mean(CY15_Table_TFDEGsintersection$CY15_24h_normalize_degree),
                                                  mean(CY16_Table_TFDEGsintersection$CY16_1h_normalize_degree),
                                                  mean(CY16_Table_TFDEGsintersection$CY16_3h_normalize_degree),
                                                  mean(CY16_Table_TFDEGsintersection$CY16_12h_normalize_degree),
                                                  mean(CY16_Table_TFDEGsintersection$CY16_24h_normalize_degree),
                                                  mean(CY20_Table_TFDEGsintersection$CY20_1h_normalize_degree),
                                                  mean(CY20_Table_TFDEGsintersection$CY20_3h_normalize_degree),
                                                  mean(CY20_Table_TFDEGsintersection$CY20_12h_normalize_degree),
                                                  mean(CY20_Table_TFDEGsintersection$CY20_24h_normalize_degree)),
                                      timecourse = rep(c("1h", "3h", "12h", "24h"), times = 3),
                                      sort = c(1:4, 1:4, 1:4)
                                      )


g <- ggplot(CY_normalizedegree_aver,
            aes(x = reorder(timecourse, sort),
                y = average,
                group = sample
                )
            )
g <- g + geom_line(aes(colour=sample))
g <- g + geom_point(aes(colour=sample))
g <- g + xlab("")
g <- g + ylab("average of normalize degree")
g <- g + theme(axis.title=element_text(size=18))
g <- g + theme(legend.text =  element_text(size = 18))

plot(g)
ggsave(file = "Nakano_RNAseq/network_analysis/results_Fig/network_degree/CY151620 TFGEGsintersection average of normalize degree.png", plot = g)
