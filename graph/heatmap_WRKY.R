T_data <- CY15_Table_3hTFDEGssort[CY15_Table_3hTFDEGssort$TF_family. == "WRKY", ][1:3, ]
data <- data.frame(normalize_degree = c(T_data$CY15_1h_normalize_degree, T_data$CY15_3h_normalize_degree, T_data$CY15_12h_normalize_degree, T_data$CY15_24h_normalize_degree),
                   time = rep(c("01h", "03h", "12h", "24h"), each = 3),
                   sample = rep(c("WRKY33", "WRKY30", "WRKY25"), times = 4)
                   )


library(ggplot2)

g <- ggplot(data, aes(x = time, y = sample, fill = normalize_degree))
g <- g + geom_tile()
g <- g + theme_bw()
g <- g + scale_fill_gradient2(low = "blue", high = "red")
g <- g + ggtitle("CY15_tGRN pathogen enrichment score")
plot(g)
ggsave("~/Nakano_RNAseq/network_analysis/results_Fig/tGRN_analysis/WRKY_relatuve_degree_heatmap.png")