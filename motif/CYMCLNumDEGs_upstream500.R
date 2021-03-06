#"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/CYMCLNumDEGs_upstream500.R"
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
t <- proc.time()
library("stringr")
CY_Table <- list(CY15_Table_20180801, CY16_Table_20180801, CY20_Table_20180801)
CY <- c("CY15", "CY16", "CY20")
time <- c("1h", "3h", "12h", "24h")
total_m <- length(CY)
m <- 1
for(m in m:total_m){
  total_n <- length(time)
  n <- 1
  for(n in n:total_n){
    temp <- CY_Table[[m]][, paste0(CY[m], "_", time[n], "_", "MCLNum")][!is.na(CY_Table[[m]][, paste0(CY[m], "_", time[n], "_", "MCLNum")])]
    T_TRUE <- c()
    T_AGI <- c()
    T_data <- c()
    T_control <- c()
    T_control_data <- c()
    total_i <- length(unique(sort(temp[temp != 0])))
    i <- 1
    for(i in i:total_i){
      T_TRUE <- intersect(which(CY_Table[[m]][, paste0(CY[m], "_", time[n], "_", "DEGs")] == "Yes"), which(i == CY_Table[[m]][, paste0(CY[m], "_", time[n], "_", "MCLNum")]))
      T_AGI <- CY_Table[[m]][T_TRUE, "AGI"]
      T_data <- up_500bp[match(T_AGI, up_500bp$AGI), ]
      T_control <- sample(up_500bp$AGI, length(T_AGI))
      T_control_data <- up_500bp[match(T_control, up_500bp$AGI), ]
      
      fasta <- c()
      allfasta <- c()
      cont_fasta <- c()
      cont_allfasta <- c()
      total_o <- nrow(T_data)
      o <- 1
      for(o in o:total_o){
        fasta <- rbind(str_sub(T_data[, "sequence"][o], start=1, end=80),
                       str_sub(T_data[, "sequence"][o], start=81, end=160),
                       str_sub(T_data[, "sequence"][o], start=161, end=240),
                       str_sub(T_data[, "sequence"][o], start=241, end=320),
                       str_sub(T_data[, "sequence"][o], start=321, end=400),
                       str_sub(T_data[, "sequence"][o], start=401, end=480),
                       str_sub(T_data[, "sequence"][o], start=481, end=500)
                       )
        cont_fasta <- rbind(str_sub(T_control_data[, "sequence"][o], start=1, end=80),
                            str_sub(T_control_data[, "sequence"][o], start=81, end=160),
                            str_sub(T_control_data[, "sequence"][o], start=161, end=240),
                            str_sub(T_control_data[, "sequence"][o], start=241, end=320),
                            str_sub(T_control_data[, "sequence"][o], start=321, end=400),
                            str_sub(T_control_data[, "sequence"][o], start=401, end=480),
                            str_sub(T_control_data[, "sequence"][o], start=481, end=500)
                            )
        data_fastaAGI <- paste0(">", T_data[, "AGI"][o])
        control_fastaAGI <- paste0(">", T_control_data[, "AGI"][o])
        allfasta <- c(allfasta, rbind(data_fastaAGI, fasta))
        cont_allfasta <- c(cont_allfasta, rbind(control_fastaAGI, cont_fasta))
        o <- o+1
      }
      target <- paste0("~/bigdata/yasue/motif/", CY[m], "/multi-fasta/MCLNum/target/", time[n], "/", CY[m], "_", time[n], "_MCLNum", i, "_upstream500.fasta")
      control <- paste0("~/bigdata/yasue/motif/", CY[m], "/multi-fasta/MCLNum/control/", time[n], "/control_", CY[m], "_", time[n], "_MCLNum", i, "_upstream500.fasta")
      write.table(allfasta, file = target, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
      write.table(cont_allfasta, file = control, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
      print(i)
      i <- i+1
    }
    n <- n+1
  }
  m <- m+1
}
t1 <- proc.time() - t
print(t1)