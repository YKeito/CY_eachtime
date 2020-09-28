sample <- c("PstDC3000_hrp.", "Botrytis.cinerea")
n <- 1
average <- c()
temp <- c("CY15_1h", "CY15_3h", "CY15_12h", "CY15_24h")
temp2 <- c("CY15_1h_q_value", "CY15_3h_q_value", "CY15_12h_q_value", "CY15_24h_q_value")
for(n in n:length(sample)){
  data <- na.omit(allRNASeq[match(CY15_Table_20180801$AGI[CY15_Table_20180801[, sample[n]] == "Yes"], rownames(allRNASeq)), ])
  i <- 1
  for(i in i:length(temp)){
    average <- c(average, mean(data[data[, temp2[i]] < 0.05, temp[i]]))
    i <- i+1
  }
  n <- n+1
}

library(ggplot2)
library(reshape2)
library(ggsci)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()

data <- data.frame(average = average,
                   sample = rep(c("PstDC3000", "B.cinerea"), each = 4),
                   time = rep(c("01h", "03h", "12h", "24h"), times = 2)
                   )


g <- ggplot(data,
            aes(x = time,
                y = average,
                group = sample
            )
)
g <- g + geom_line(aes(colour=sample))
g <- g + geom_point(aes(colour=sample))
g <- g + xlab("")
g <- g + ylab("Average expression")
g <- g+ylim(c(0, NA))

plot(g)