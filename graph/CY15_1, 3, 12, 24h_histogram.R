object <- list(CY15_Table_1hTFDEGssort, CY15_Table_3hTFDEGssort, CY15_Table_12hTFDEGssort, CY15_Table_24hTFDEGssort, CY15_Table_TFDEGsintersection)
Figtitle <- c("CY15_Table_1h_TFDEGs", "CY15_Table_3h_TFDEGs", "CY15_Table_12h_TFDEGs", "CY15_Table_24h_TFDEGs", "CY15_Table_TFDEGs_intersection")
g <-c()
library(ggplot2)
#パッケージのインストール
#install.packages("devtools")
#devtools::install_github("dgrtwo/gganimate")
library("gganimate")

i <- 1
for(i in i:length(object)){
  AGI <- rep(object[[i]][, "AGI"], times = 4)
  timecourse <- rep(c("01h", "03h", "12h", "24h"), each = nrow(object[[i]]))
  normalize_degree <- c(object[[i]][, "CY15_1h_normalize_degree"],
                        object[[i]][, "CY15_3h_normalize_degree"],
                        object[[i]][, "CY15_12h_normalize_degree"],
                        object[[i]][, "CY15_24h_normalize_degree"]
  )
  data <- data.frame(AGI = AGI,
                     timecourse = timecourse,
                     normalize_degree = normalize_degree
                     )
  
  levels(data$timecourse) <- c("1h", "3h", "12h", "24h")
  g <- ggplot(data, aes(x = normalize_degree, fill = timecourse))
  g <- g + geom_histogram(position = "identity", alpha = 0.6, binwidth=0.1)
  g <- g + facet_wrap(timecourse ~ ., nrow=5) 
  g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
  Figtitle1 <- paste0(Figtitle[i], "NodeNum:", nrow(object[[i]]), " genes")
  g <- g + ggtitle(Figtitle1)
  plot(g)
  #plot(g)
  output <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/network_degree/histogram/TFDEGs/", Figtitle1, "_histogram", ".png")
  ggsave(file = output, plot = g)
  i <- i+1
}
