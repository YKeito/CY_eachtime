MCLClu_1h <- unique(CY15_Table_20180801$CY15_1h_MCLNum)
MCLClu_1h <- MCLClu_1h[order(MCLClu_1h, decreasing = F)]
temp <- c()
i <- 1
for(i in i:length(MCLClu_1h)){
  temp <- c(temp, sum(CY15_Table_20180801$CY15_1h_MCLNum == MCLClu_1h[i], na.rm = T))
}
temp <- temp[2:length(temp)]


MCLClu_3h <- unique(CY15_Table_20180801$CY15_3h_MCLNum)
MCLClu_3h <- MCLClu_3h[order(MCLClu_3h, decreasing = F)]
temp2 <- c()
i <- 1
for(i in i:length(MCLClu_3h)){
  temp2 <- c(temp2, sum(CY15_Table_20180801$CY15_3h_MCLNum == MCLClu_3h[i], na.rm = T))
}
temp2 <- temp2[2:length(temp2)]


MCLClu_12h <- unique(CY15_Table_20180801$CY15_12h_MCLNum)
MCLClu_12h <- MCLClu_12h[order(MCLClu_12h, decreasing = F)]

temp3 <- c()
i <- 1
for(i in i:length(MCLClu_12h)){
  temp3 <- c(temp3, sum(CY15_Table_20180801$CY15_12h_MCLNum == MCLClu_12h[i], na.rm = T))
}
temp3 <- temp3[2:length(temp3)]

MCLClu_24h <- unique(CY15_Table_20180801$CY15_24h_MCLNum)
MCLClu_24h <- MCLClu_24h[order(MCLClu_24h, decreasing = F)]
temp4 <- c()
i <- 1
for(i in i:length(MCLClu_24h)){
  temp4 <- c(temp4, sum(CY15_Table_20180801$CY15_24h_MCLNum == MCLClu_24h[i], na.rm = T))
}
temp4 <- temp4[2:length(temp4)]


len <- c(length(temp), length(temp2), length(temp3), length(temp4))
data <- matrix(rep(0, times = max(len)), ncol = 4, nrow = max(len))
n <- 1
total <- max(len)
T_all <- c()
for(n in n:total){
  Num <- "000"
  Num <- paste0(substr(Num, 1, nchar(Num)-nchar(n)), n)
  T_all <- c(T_all, Num)
  print(n)
  n <- n+1
}
colnames(data) <- c("1h", "3h", "12h", "24h")
rownames(data) <- T_all
data[1:length(temp), "1h"] <- temp
data[1:length(temp2), "3h"] <- temp2
data[1:length(temp3), "12h"] <- temp3
data[1:length(temp4), "24h"] <- temp4

g <- data.frame(total = c(data[, 1], data[, 2], data[, 3], data[, 4]),
                NumMCL = rep(T_all, times = 4),
                time = rep(c("01h", "03h", "12h", "24h"), each = nrow(data))
                )
library(ggplot2)
g <- ggplot(g,
            aes(x = NumMCL,
                y = total,
                group = time
            )
)
g <- g + geom_line(aes(colour=time))
g <- g + geom_point(aes(colour=time))

plot(g)
ggsave(file = "Nakano_RNAseq/network_analysis/results_Fig/network_degree/CY151620 TFGEGs average of normalize degree.png", plot = g)
