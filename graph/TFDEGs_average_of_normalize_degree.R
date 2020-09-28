library(ggplot2)
library(reshape2)
CY_normalizedegree_aver <- data.frame(sample = rep(c("CY15", "CY16", "CY20"), each = 4),
                                      average = c(mean(CY15_Table_1hTFDEGssort$CY15_1h_normalize_degree[CY15_Table_1hTFDEGssort$CY15_1h_normalize_degree != 0]), 
                                                  mean(CY15_Table_3hTFDEGssort$CY15_3h_normalize_degree[CY15_Table_3hTFDEGssort$CY15_3h_normalize_degree != 0]),
                                                  mean(CY15_Table_12hTFDEGssort$CY15_12h_normalize_degree[CY15_Table_12hTFDEGssort$CY15_12h_normalize_degree != 0]),
                                                  mean(CY15_Table_24hTFDEGssort$CY15_24h_normalize_degree[CY15_Table_24hTFDEGssort$CY15_24h_normalize_degree != 0]),
                                                  mean(CY16_Table_1hTFDEGssort$CY16_1h_normalize_degree[CY16_Table_1hTFDEGssort$CY16_1h_normalize_degree != 0]),
                                                  mean(CY16_Table_3hTFDEGssort$CY16_3h_normalize_degree[CY16_Table_3hTFDEGssort$CY16_3h_normalize_degree != 0]),
                                                  mean(CY16_Table_12hTFDEGssort$CY16_12h_normalize_degree[CY16_Table_12hTFDEGssort$CY16_12h_normalize_degree != 0]),
                                                  mean(CY16_Table_24hTFDEGssort$CY16_24h_normalize_degree[CY16_Table_24hTFDEGssort$CY16_24h_normalize_degree != 0]),
                                                  mean(CY20_Table_1hTFDEGssort$CY20_1h_normalize_degree[CY20_Table_1hTFDEGssort$CY20_1h_normalize_degree != 0]),
                                                  mean(CY20_Table_3hTFDEGssort$CY20_3h_normalize_degree[CY20_Table_3hTFDEGssort$CY20_3h_normalize_degree != 0]),
                                                  mean(CY20_Table_12hTFDEGssort$CY20_12h_normalize_degree[CY20_Table_12hTFDEGssort$CY20_12h_normalize_degree != 0]),
                                                  mean(CY20_Table_24hTFDEGssort$CY20_24h_normalize_degree[CY20_Table_24hTFDEGssort$CY20_24h_normalize_degree != 0])
                                                  )
                                      ,
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
ggsave(file = "Nakano_RNAseq/network_analysis/results_Fig/network_degree/CY151620 TFGEGs average of normalize degree.png", plot = g)



