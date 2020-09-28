#"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/motif/CYMCLNumDEGs_motif_pvalue.R"
#library("stringr")
#allmotif <- read.table("~/bigdata/yasue/MEME/Databases/motif_databases/CIS-BP/Arabidopsis_thaliana.meme", sep = "\t", header = T, fill = T, stringsAsFactors = F)
all <- allmotif$MEME.version.4[grep("MOTIF", allmotif$MEME.version.4)]
all <- str_sub(all, start = 7, end = 16)

time <- c("1h", "3h", "12h", "24h")
T_pvalue <- c()
obname <- c()
m <- 1
for(m in m:length(allMotif_ID)){
  temp <- allMotif_ID[[m]][!is.na(allMotif_ID[[m]])]
  target <- unique(names(temp))
  n <- 1
  for(n in n:length(target)){
    i <- 1
    for(i in i:length(unique(all))){
      nn <- length(temp[names(temp) == target[n]])
      a <- length(grep(unique(all)[i], temp[names(temp) == target[n]], invert = F))
      MM <- sum(all == unique(all)[i])
      b <- MM-a
      c <- nn-a
      d <- length(all)-a-b-c
      NN <-length(all)
      mx <- matrix(c(a, b, c, d), nrow=2, byrow=T)
      T_pvalue <- c(T_pvalue, fisher.test(mx, alternative = "greater")$p.value)
      obname <- c(obname, paste0("CY15_", time[m], "_", target[n], "_MotifID:", all[i]))
      i <- i+1
    }
    n <- n+1
  }
  print(m)
  m <- m+1
}
names(T_pvalue) <- obname
length(T_pvalue[which(T_pvalue < 0.05)])
