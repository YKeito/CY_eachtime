temp <- list(c(CY20_Table_20180801$TF_family. != "No" & CY20_Table_20180801$CY20_1h_DEGs == "Yes" & CY20_Table_20180801$CY20_1h_normalize_degree != 0 & CY20_Table_20180801$CY20_3h_normalize_degree != 0 & CY20_Table_20180801$CY20_12h_normalize_degree != 0 & CY20_Table_20180801$CY20_24h_normalize_degree != 0),
             c(CY20_Table_20180801$TF_family. != "No" & CY20_Table_20180801$CY20_3h_DEGs == "Yes" & CY20_Table_20180801$CY20_1h_normalize_degree != 0 & CY20_Table_20180801$CY20_3h_normalize_degree != 0 & CY20_Table_20180801$CY20_12h_normalize_degree != 0 & CY20_Table_20180801$CY20_24h_normalize_degree != 0),
             c(CY20_Table_20180801$TF_family. != "No" & CY20_Table_20180801$CY20_12h_DEGs == "Yes" & CY20_Table_20180801$CY20_1h_normalize_degree != 0 & CY20_Table_20180801$CY20_3h_normalize_degree != 0 & CY20_Table_20180801$CY20_12h_normalize_degree != 0 & CY20_Table_20180801$CY20_24h_normalize_degree != 0),
             c(CY20_Table_20180801$TF_family. != "No" & CY20_Table_20180801$CY20_24h_DEGs == "Yes" & CY20_Table_20180801$CY20_1h_normalize_degree != 0 & CY20_Table_20180801$CY20_3h_normalize_degree != 0 & CY20_Table_20180801$CY20_12h_normalize_degree != 0 & CY20_Table_20180801$CY20_24h_normalize_degree != 0)
)
temp2 <- c("AGI", "TF_family.", colnames(CY20_Table_20180801)[grep("normalize_degree", colnames(CY20_Table_20180801))])
temp3 <- c("CY20_1h_normalize_degree", "CY20_3h_normalize_degree", "CY20_12h_normalize_degree", "CY20_24h_normalize_degree")
title <- c("CY20 1h", "CY20 3h", "CY20 12h", "CY20 24h")

n <- 1
for(n in n:length(temp)){
  data <- CY20_Table_20180801[temp[[n]], temp2]
  data <- rbind(c(NA, NA, rep(1, times = 4)), data, c(NA, NA, rep(0, times = 4)))
  
  AGI <- rep(data[, "AGI"], times = 4)
  timecourse <- rep(c("01h", "03h", "12h", "24h"), each = nrow(data))
  sort <- order(data[, temp3[n]], decreasing = F)
  normalize_degree <- c(data[sort, temp3[1]],
                        data[sort, temp3[2]],
                        data[sort, temp3[3]],
                        data[sort, temp3[4]]
  )
  
  gg <- data.frame(AGI = AGI,
                   timecourse = timecourse,
                   normalize_degree = normalize_degree,
                   reorder = rep(1:length(sort), times = 4),
                   reordertime = rep(c(1:4), each = nrow(data)),
                   stringsAsFactors = F
  )
  g <- ggplot(gg, aes(x = reorder(timecourse, reordertime), y = reorder(AGI, reorder), fill = normalize_degree))
  g <- g + geom_tile()
  g <- g + theme_bw()
  g <- g + scale_fill_gradient2(low = "blue", high = "red", na.value = "black")
  
  Figtitle1 <- paste0(title[n], " NodeNum:", nrow(data), " TFDEGs")
  g <- g + ggtitle(Figtitle1)
  output <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/network_degree/CY20/heatmap/all_TFDEGs/", Figtitle1, " heatmap", ".png")
  ggsave(file = output, plot = g, dpi = 100, width = 8, height = 6)
  plot(g)
  print(n)
  n <- n+1
}