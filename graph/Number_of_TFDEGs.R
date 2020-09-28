library(ggplot2)
library(reshape2)
CY_normalizedegree_aver <- data.frame(sample = rep(c("CY15", "CY16", "CY20"), each = 4),
                                      average = c(nrow(CY15_Table_1hTFDEGssort), 
                                                  nrow(CY15_Table_3hTFDEGssort),
                                                  nrow(CY15_Table_12hTFDEGssort),
                                                  nrow(CY15_Table_24hTFDEGssort),
                                                  nrow(CY16_Table_1hTFDEGssort),
                                                  nrow(CY16_Table_3hTFDEGssort),
                                                  nrow(CY16_Table_12hTFDEGssort),
                                                  nrow(CY16_Table_24hTFDEGssort),
                                                  nrow(CY20_Table_1hTFDEGssort),
                                                  nrow(CY20_Table_3hTFDEGssort),
                                                  nrow(CY20_Table_12hTFDEGssort),
                                                  nrow(CY20_Table_24hTFDEGssort)
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
g <- g + ylab("Number of TFs")
g <- g + theme(axis.title=element_text(size=18))
g <- g + theme(legend.text =  element_text(size = 18))

plot(g)
ggsave(file = "Nakano_RNAseq/network_analysis/results_Fig/network_degree/Number of TF:CY151620 TFGEGs.png", plot = g)



