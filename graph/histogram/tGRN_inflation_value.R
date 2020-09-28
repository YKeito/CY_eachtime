#"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/graph/histogram/tGRN_inflation_value.R"
filename <- list.files("~/Nakano_RNAseq/network_analysis/base/CY151620_131224h_mastercluster_cytoscapeformat", pattern = ".csv", full.names = T)
object_name <- c("inflation2.5", "inflation3", "inflation4", "inflation5")
temp <- list(inflation2.5, inflation3, inflation4, inflation5)
MCLNum <- c()
Category <- c()
library(ggplot2)
n <- 1
for (n in 1:length(filename)){
  assign(object_name[n], read.table(filename[n],  header=T, sep=",", stringsAsFactors = F))
  MCLNum <- c(temp[[n]][, "X__mclCluster"][temp[[n]][, "X__mclCluster"] < 6][!is.na(temp[[n]][, "X__mclCluster"][temp[[n]][, "X__mclCluster"] < 6])])
  Category <- rep(object_name[n], times = length(MCLNum[[n]]))
  T_data <- data.frame(MCLNum = unlist(MCLNum),
                       Category = Category
  )
  g <- ggplot(T_data, aes(x = MCLNum))
  g <- g + geom_histogram(binwidth=0.5)
  g <- g + ggtitle(object_name[n])
  plot(g)
  title <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/tGRN_inflation_histogram/", object_name[n], ".png")
  ggsave(title, g)
  n <- n+1
}
