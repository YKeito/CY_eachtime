#"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/CYDEGs_upstream500.R"
###########################################
#upstream_500 <- read.table(file = "~/bigdata/yasue/motif/TAIR10_upstream_500_20101028", fill = T, sep = ",", stringsAsFactors = F)
#upstream_500 <- as.character(unlist(upstream_500))
temp <- grep("chr", upstream_500)

AGI <- c()
i <- 1
total <- length(temp)
for(i in i:total){
  AGI <- c(AGI, substr(upstream_500[temp[i]], 2, 10))
  print(i)
  i <- i+1
}

#最後の一周だけ自動化できなかった。
allsequence <- list()
presequence <- c()
sequence <- c()
i <- 1
total <- length(temp)
for(i in i:c(total-1)){
  presequence <- upstream_500[c(temp[i]+1):c(temp[i+1]-1)]
  n <- 1
  for(n in n:length(presequence)){
    sequence <- paste0(sequence, presequence[n])
    n <- n+1
  }
  allsequence <- c(allsequence, list(sequence))
  sequence <- c()
  print(i)
  i <- i+1
}
#残りの最後の一周を追加
presequence <- upstream_500[c(temp[total]+1):length(upstream_500)]
sequence <- paste0(sequence, presequence[n])
allsequence <- c(allsequence, list(sequence))

#data.frame
up_500bp <- data.frame(AGI = AGI,
                       sequence = unlist(allsequence),
                       stringsAsFactors = F
)
#####################対象の遺伝子群を引っ張ってくる#######################
temp <- CY20_Table_20180801[, "CY20_24h_DEGs"] == "Yes"
temp <- CY20_Table_20180801[, "AGI"][temp]
temp <- na.omit(temp)
CY20_sequence <- up_500bp[match(temp, up_500bp$AGI), ]
control <- sample(up_500bp$AGI, length(temp))
control_sequence <- up_500bp[match(control, up_500bp$AGI), ]
#txt -> fasta
library("stringr")
i <- 1
fasta <- c()
allfasta <- c()

cont_fasta <- c()
cont_allfasta <- c()


total <- nrow(CY20_sequence)
for(i in i:total){
  fasta <- rbind(str_sub(CY20_sequence$sequence[i], start=1, end=80),
                 str_sub(CY20_sequence$sequence[i], start=81, end=160),
                 str_sub(CY20_sequence$sequence[i], start=161, end=240),
                 str_sub(CY20_sequence$sequence[i], start=241, end=320),
                 str_sub(CY20_sequence$sequence[i], start=321, end=400),
                 str_sub(CY20_sequence$sequence[i], start=401, end=480),
                 str_sub(CY20_sequence$sequence[i], start=481, end=500)
  )
  cont_fasta <- rbind(str_sub(control_sequence$sequence[i], start=1, end=80),
                      str_sub(control_sequence$sequence[i], start=81, end=160),
                      str_sub(control_sequence$sequence[i], start=161, end=240),
                      str_sub(control_sequence$sequence[i], start=241, end=320),
                      str_sub(control_sequence$sequence[i], start=321, end=400),
                      str_sub(control_sequence$sequence[i], start=401, end=480),
                      str_sub(control_sequence$sequence[i], start=481, end=500)
                      )
  
  CY20_fastaAGI <- paste0(">", CY20_sequence$AGI[i])
  control_fastaAGI <- paste0(">", control_sequence$AGI[i])
  
  allfasta <- c(allfasta, rbind(CY20_fastaAGI, fasta))
  cont_allfasta <- c(cont_allfasta, rbind(control_fastaAGI, cont_fasta))
  
  print(i)
  i <- i+1
}

write.table(allfasta, "~/bigdata/yasue/motif/CY20/multi-fasta/CY20_24h_upstream500.fasta", append = F, quote = F, sep = "\t", row.names = F, col.names = F)
write.table(cont_allfasta, "~/bigdata/yasue/motif/CY20/multi-fasta/control_CY20_24h_upstream500.fasta", append = F, quote = F, sep = "\t", row.names = F, col.names = F)