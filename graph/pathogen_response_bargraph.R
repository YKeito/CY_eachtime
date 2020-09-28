#"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/graph/heatmap_pathogen_response.R"
library(ggplot2)

Botrytis_cinerea <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/B.cinerea_newlist.txt", sep = "\t", header = T, stringsAsFactors = F)
Botrytis_cinerea$Gene.Locus <- toupper(Botrytis_cinerea$Gene.Locus)
test <- c("PstDC3000_hrp.", "Botrytis.cinerea")
sample <- list(PstDC3000_hrp$AGI, Botrytis_cinerea$Gene.Locus)
CY <- list(CY15_Table_20180801, CY16_Table_20180801, CY20_Table_20180801)
temp <- c("CY15", "CY16", "CY20")
T_pvalue <- c()
obnames <- c()
n <- 1
for(n in n:length(temp)){
  time <- paste0(temp[n], c("_1h_DEGs", "_3h_DEGs", "_12h_DEGs", "_24h_DEGs"))
  m <- 1
  for(m in m:length(sample)){
    i <- 1
    for(i in i:length(time)){
      a <- length(intersect(sample[[m]], CY[[n]][, "AGI"][CY[[n]][, time[[i]]] == "Yes"]))
      b <- length(setdiff(sample[[m]], CY[[n]][, "AGI"][CY[[n]][, time[[i]]] == "Yes"]))
      c <- length(setdiff(CY[[n]][, "AGI"][CY[[n]][, time[[i]]] == "Yes"], sample[[m]]))
      d <- sum(CY[[n]][, test[m]] == "No" & CY[[n]][, time[[i]]] == "No")
      mx <- matrix(c(a, b, c, d), nrow=2, byrow=T)
      T_pvalue <- c(T_pvalue, fisher.test(mx, alternative = "greater")$p.value)
      i <- i+1
    }
    obnames <- c(obnames, paste0(time, ",",test[m]))
    m <- m+1
  }
}

data <- data.frame(sample = rep(rep(test, each = 4), times = 3),
                   timecourse = rep(c("01h", "03h", "12h", "24h"), times = 2*3),
                   CY = rep(temp, each = 8),
                   enrichment_score = -log10(T_pvalue)
)

g <- ggplot(data, aes(x = timecourse, y = enrichment_score, fill = sample))
g <- g + geom_bar(stat = "identity", position = "dodge")
g <- g +  theme_bw()
g <- g + facet_wrap(CY ~ ., ncol=5)
g <- g + scale_fill_grey(start = .2, end = .7, guide = "none")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
g <- g + theme(legend.position = 'none')
plot(g)
title <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/tGRN_analysis/pathogen_enrichment_Fig/", "CYall", "_tGRN_pathogen_enrichment_bargraph.png")
ggsave(title, width = 8, height = 4)