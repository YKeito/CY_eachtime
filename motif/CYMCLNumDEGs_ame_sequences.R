#*tGRN_MCL:object:"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/graph/tGRN_MCL_AverageOfExpression.R"
#this script:"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/motif/CYMCLNumDEGs_ame_sequences.R"
temp <- CY15_tGRN_MCL[CY15_tGRN_MCL$MCLClu != 2 & CY15_tGRN_MCL$MCLClu != 1, ]

time <- c("1h", "3h", "12h", "24h")
allMotif_ID <- list()
n <- 1
for(n in n:length(time)){
  T_data <- temp[grep(time[n], rownames(temp)), ]
  title <- paste0("~/bigdata/yasue/motif/CY15/motif_results/", time[n])
  T_filename <- list.files(title, full.names = T)
  T_MCLNum <- as.numeric(T_data$MCLNum)
  filename <- c()
  Motif_ID <- c()
  obnames <- c()
  i <- 1
  for(i in i:length(T_MCLNum)){
    filename <- T_filename[grep(paste0("MCLNum", T_MCLNum[i], "_"), T_filename)]
    title <- paste0(filename, "/sequences.tsv")
    ame_sequences_results <- try(read.table(title, sep = "\t", header = T, stringsAsFactors = F), silent = TRUE)
    if(class(ame_sequences_results) != "try-error"){
      Motif_ID <- c(Motif_ID, ame_sequences_results$motif_ID)
      obnames <- c(obnames, rep(paste0("MCLNum", T_MCLNum[i]), times = length(ame_sequences_results$motif_ID)))
    }else{
      Motif_ID <- c(Motif_ID, NA)
      obnames <- c(obnames, paste0("MCLNum", T_MCLNum[i]))
    }
    i <- i+1
  }
  names(Motif_ID) <- obnames
  allMotif_ID <- c(allMotif_ID, list(Motif_ID))
  print(n)
  n <- n+1
}
