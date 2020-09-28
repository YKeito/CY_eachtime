#*tGRN_MCL:object:"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/graph/tGRN_MCL_AverageOfExpression.R"
#this script:"~/Nakano_RNAseq/network_analysis/script/CY_eachtime/motif/CYMCLNumDEGs_motif.R"
#motif <- read.table(file = "~/bigdata/yasue/motif/AtCisDB/BindingSite.tbl", sep = "\t", header = T, fill = T, stringsAsFactors = F)
temp <- CY15_tGRN_MCL[CY15_tGRN_MCL$MCLClu != 2 & CY15_tGRN_MCL$MCLClu != 1, ]

time <- c("1h", "3h", "12h", "24h")
allMotif_ID <- list()
allconsensus <- list()
n <- 1
for(n in n:length(time)){
  T_data <- temp[grep(time[n], rownames(temp)), ]
  title <- paste0("~/bigdata/yasue/motif/CY15/motif_results/", time[n])
  T_filename <- list.files(title, full.names = T)
  T_MCLNum <- as.numeric(T_data$MCLNum)
  filename <- c()
  Motif_ID <- c()
  obnames <- c()
  consensus <- c()
  i <- 1
  for(i in i:length(T_MCLNum)){
    filename <- T_filename[grep(paste0("MCLNum", T_MCLNum[i], "_"), T_filename)]
    title <- paste0(filename, "/ame.tsv")
    motif_ame_results <- try(read.table(title, sep = "\t", header = T, stringsAsFactors = F), silent = TRUE)
    if(class(motif_ame_results) != "try-error"){
      Motif_ID <- c(Motif_ID, motif_ame_results$motif_ID)
      consensus <- c(consensus, motif_ame_results$consensus)
      obnames <- c(obnames, rep(paste0("MCLNum", T_MCLNum[i]), times = length(motif_ame_results$motif_ID)))
    }else{
      Motif_ID <- c(Motif_ID, NA)
      consensus <- c(consensus, NA)
      obnames <- c(obnames, paste0("MCLNum", T_MCLNum[i]))
    }
    i <- i+1
  }
  names(Motif_ID) <- obnames
  allMotif_ID <- c(allMotif_ID, list(Motif_ID))
  allconsensus <- c(allconsensus, list(consensus))
  print(n)
  n <- n+1
}

CY15_MotifID_consensus <- data.frame(Name = c(paste0("CY15_1h_", names(unlist(allMotif_ID[1]))), 
                                              paste0("CY15_3h_", names(unlist(allMotif_ID[2]))), 
                                              paste0("CY15_12h_", names(unlist(allMotif_ID[3]))), 
                                              paste0("CY15_24h_", names(unlist(allMotif_ID[4])))),
                                     MotifID = unlist(allMotif_ID),
                                     consensu = unlist(allconsensus),
                                     stringsAsFactors = F
)

CY15_MotifID_consensus <- na.omit(CY15_MotifID_consensus)

data <- data.frame(sample = rep(unique(CY15_MotifID_consensus$Name), each = length(unique(CY15_MotifID_consensus$MotifID))),
                   MotifID = rep(unique(CY15_MotifID_consensus$MotifID), times = length(unique(CY15_MotifID_consensus$Name))),
                   value = rep(0, times = length(unique(CY15_MotifID_consensus$Name))*length(unique(CY15_MotifID_consensus$MotifID))),
                   stringsAsFactors = F
                   )

target <- c(paste0("CY15_1h_", names(unlist(allMotif_ID[[1]])[!is.na(unlist(allMotif_ID[[1]]))]), unlist(allMotif_ID[[1]])[!is.na(unlist(allMotif_ID[[1]]))]),
             paste0("CY15_3h_", names(unlist(allMotif_ID[[2]])[!is.na(unlist(allMotif_ID[[2]]))]), unlist(allMotif_ID[[2]])[!is.na(unlist(allMotif_ID[[2]]))]),
             paste0("CY15_12h_", names(unlist(allMotif_ID[[3]])[!is.na(unlist(allMotif_ID[[3]]))]), unlist(allMotif_ID[[3]])[!is.na(unlist(allMotif_ID[[3]]))]),
             paste0("CY15_24h_", names(unlist(allMotif_ID[[4]])[!is.na(unlist(allMotif_ID[[4]]))]), unlist(allMotif_ID[[4]])[!is.na(unlist(allMotif_ID[[4]]))])
             )
temp <- paste0(data$sample, data$MotifID)
i <- 1
for(i in i:length(target)){
  data$value[which(target[i] == temp)] <- 1
}

sample_sort <- c()
n <- 1
for(n in n:length(unique(data$sample))){
  sample_sort <- c(sample_sort, rep(c(length(unique(data$sample))+1-n), times = sum(data$sample == unique(data$sample)[n])))
}

df <- data.frame(data, 
                 sample_sort,
                 stringsAsFactors = F
                 )

g <- ggplot(df, aes(x = MotifID, y = reorder(sample, sample_sort), fill = value))
g <- g + geom_tile(color = "black", size = 0.1)
g <- g + scale_fill_gradient2(high = "red")
g <- g+theme_dark()
g <- g + theme_linedraw()
#g <- g + coord_flip()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(g)
ggsave("~/Nakano_RNAseq/network_analysis/results_Fig/tGRN_Subcluster_MotifID_Enrichment/CY15/CY15_tGRN_Subcluster_MotifID_AME_heatmap.png", g, width = 20, height = 10)