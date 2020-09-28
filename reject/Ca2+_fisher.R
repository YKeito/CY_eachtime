Calcium <- read.table(file = "~/Nakano_RNAseq/network_analysis/base/genes_set/calcium-responsive genes/calcium-responsive genes.txt", sep = "\t", stringsAsFactors = F, header = T)
total <- nrow(allRNASeq_foldchange)
intersection <- length(intersect(Calcium$AGI, rownames(CY15_1h_DEGs)))
a <- c(nrow(CY15_1h_DEGs)-intersection)
b <- nrow(Calcium)-intersection
c <- c(total-nrow(Calcium)-a)

test <- matrix(c(intersection, a, b, c), ncol = 2, nrow = 2)
fisher.test(test)$p.value


x <- length(intersect(Calcium$AGI, rownames(CY15_1h_DEGs)))
M <- nrow(Calcium)
N <- nrow(allRNASeq_foldchange)
n <- nrow(CY15_1h_DEGs)

phyper(c(x-1), M, c(N-M), n, lower.tail = F)

