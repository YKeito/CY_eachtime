library(ggplot2)
object <- list(CY15_Table_1hTFDEGssort, CY15_Table_3hTFDEGssort, CY15_Table_12hTFDEGssort, CY15_Table_24hTFDEGssort,
               CY15_Table_TFDEGsintersection)
Figtitle <- c("CY15_Table_1h_TFDEGs", "CY15_Table_3h_TFDEGs", "CY15_Table_12h_TFDEGs", "CY15_Table_24h_TFDEGs", "CY15_Table_TFDEGs_intersection")
g <-c()
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
  g1 <- ggplot(data, aes(x = timecourse, y = normalize_degree, fill = timecourse))
  g1 <- g1+geom_boxplot(alpha=0.5,colour="gray30")
  Figtitle1 <- paste0(Figtitle[i], "NodeNum:", nrow(object[[i]]))
  g1 <- g1 + ggtitle(Figtitle1)
  plot(g1)
  output <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/network_degree/boxplot/TFDEGs/", Figtitle1, ".png")
  ggsave(file = output, plot = g1)
  i <- i+1
}
