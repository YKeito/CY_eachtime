####FDR005####
#up
CY15_1h_FDR005_up <- sum(allRNASeq$CY15_1h[allRNASeq$CY15_1h_q_value < 0.05] > 0)
CY15_3h_FDR005_up <- sum(allRNASeq$CY15_3h[allRNASeq$CY15_3h_q_value < 0.05] > 0)
CY15_12h_FDR005_up <- sum(allRNASeq$CY15_12h[allRNASeq$CY15_12h_q_value < 0.05] > 0)
CY15_24h_FDR005_up <- sum(allRNASeq$CY15_24h[allRNASeq$CY15_24h_q_value < 0.05] > 0)

CY16_1h_FDR005_up <- sum(allRNASeq$CY16_1h[allRNASeq$CY16_1h_q_value < 0.05] > 0)
CY16_3h_FDR005_up <- sum(allRNASeq$CY16_3h[allRNASeq$CY16_3h_q_value < 0.05] > 0)
CY16_12h_FDR005_up <- sum(allRNASeq$CY16_12h[allRNASeq$CY16_12h_q_value < 0.05] > 0)
CY16_24h_FDR005_up <- sum(allRNASeq$CY16_24h[allRNASeq$CY16_24h_q_value < 0.05] > 0)

CY20_1h_FDR005_up <- sum(allRNASeq$CY20_1h[allRNASeq$CY20_1h_q_value < 0.05] > 0)
CY20_3h_FDR005_up <- sum(allRNASeq$CY20_3h[allRNASeq$CY20_3h_q_value < 0.05] > 0)
CY20_12h_FDR005_up <- sum(allRNASeq$CY20_12h[allRNASeq$CY20_12h_q_value < 0.05] > 0)
CY20_24h_FDR005_up <- sum(allRNASeq$CY20_24h[allRNASeq$CY20_24h_q_value < 0.05] > 0)

#down
CY15_1h_FDR005_down <- sum(allRNASeq$CY15_1h[allRNASeq$CY15_1h_q_value < 0.05] < 0)
CY15_3h_FDR005_down <- sum(allRNASeq$CY15_3h[allRNASeq$CY15_3h_q_value < 0.05] < 0)
CY15_12h_FDR005_down <- sum(allRNASeq$CY15_12h[allRNASeq$CY15_12h_q_value < 0.05] < 0)
CY15_24h_FDR005_down <- sum(allRNASeq$CY15_24h[allRNASeq$CY15_24h_q_value < 0.05] < 0)

CY16_1h_FDR005_down <- sum(allRNASeq$CY16_1h[allRNASeq$CY16_1h_q_value < 0.05] < 0)
CY16_3h_FDR005_down <- sum(allRNASeq$CY16_3h[allRNASeq$CY16_3h_q_value < 0.05] < 0)
CY16_12h_FDR005_down <- sum(allRNASeq$CY16_12h[allRNASeq$CY16_12h_q_value < 0.05] < 0)
CY16_24h_FDR005_down <- sum(allRNASeq$CY16_24h[allRNASeq$CY16_24h_q_value < 0.05] < 0)

CY20_1h_FDR005_down <- sum(allRNASeq$CY20_1h[allRNASeq$CY20_1h_q_value < 0.05] < 0)
CY20_3h_FDR005_down <- sum(allRNASeq$CY20_3h[allRNASeq$CY20_3h_q_value < 0.05] < 0)
CY20_12h_FDR005_down <- sum(allRNASeq$CY20_12h[allRNASeq$CY20_12h_q_value < 0.05] < 0)
CY20_24h_FDR005_down <- sum(allRNASeq$CY20_24hh[allRNASeq$CY20_24h_q_value < 0.05] < 0)

#CYall FDR005
df <- data.frame(expression_change = rep(c("01up", "02down"), times = 12),
                 Numgenes = c(CY15_1h_FDR005_up, CY15_1h_FDR005_down,
                              CY15_3h_FDR005_up, CY15_3h_FDR005_down,
                              CY15_12h_FDR005_up, CY15_12h_FDR005_down,
                              CY15_24h_FDR005_up, CY15_24h_FDR005_down,
                              CY16_1h_FDR005_up, CY16_1h_FDR005_down,
                              CY16_3h_FDR005_up, CY16_3h_FDR005_down,
                              CY16_12h_FDR005_up, CY16_12h_FDR005_down,
                              CY16_24h_FDR005_up, CY16_24h_FDR005_down,
                              CY20_1h_FDR005_up, CY20_1h_FDR005_down,
                              CY20_3h_FDR005_up, CY20_3h_FDR005_down,
                              CY20_12h_FDR005_up, CY20_12h_FDR005_down,
                              CY20_24h_FDR005_up, CY20_24h_FDR005_down
                              ),
                 Category = rep(rep(c("011h", "023h", "0312h", "0424h"), each = 2), times = 3),
                 CY = rep(c("CY15", "CY16", "CY20"), each = 8)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()

levels(df$Category) <- c("1h", "3h", "12h", "24h")
levels(df$expression_change) <- c("up", "down")
g <- ggplot(
  df,
  aes (
    x = Category,             # 遺伝子別でグルーピング
    y = Numgenes,
    fill = expression_change       # 縦軸を生物種で fill up
  )
)
g <- g + geom_bar(stat = "identity")
g <- g +  theme_bw()
g <- g + scale_colour_manual(values = "Set1")
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
CYall <- g + facet_grid(CY ~ .)
CYall <- CYall + theme(legend.position = 'none')
plot(CYall)