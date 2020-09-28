#################################################################################################set##########################################################################################
#TF_family <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/arabidopsis_TF_family.txt", sep = "\t", header = T, stringsAsFactors = F)
#Botrytis_cinerea <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/pathogen infection/B.cinerea.txt", sep = "\t", header = T, stringsAsFactors = F)
#PstDC3000 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/pathogen infection/pseudomonas syringae.txt", sep = "\t", header = T, stringsAsFactors = F)
#MeJA_DEGs <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/MeJA_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
#BTH <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/BTH_affected_AGI.txt", sep = "\t", header = T, stringsAsFactors = F)
#CYall <- read.table("~/Nakano_RNAseq/network_analysis/base/CY_base/CY151620_notsignificant.txt", sep = "\t", header = T, stringsAsFactors = F)
CY_unionAGI <- CYall$gene_id

timecourse <- c("1h", "3h", "12h", "24h")
b <- c("sub1","Botrytis_cinerea", "PstDC3000", "MeJA", "BTH")
c <- c("DEGs", "degree", "BetweennessCentrality")
p <- 1
cname <- c()
for(p in p:length(timecourse)){
  n <- 1
  for(n in n:length(b)){
    o <- 1
    for(o in o:length(c)){
      cname <- c(cname, paste0("CY15_", timecourse[p], "_", b[n], "_", c[o]))
      o <- o+1
    }
    n <- n+1
  }
  p <- p+1
}

CY15_Table <- matrix(rep("No", times = length(CYall$gene_id)*62), nrow = length(CYall$gene_id), ncol = 62)
colnames(CY15_Table) <- c("AGI", "TF", cname)
colCY15_Table <- colnames(CY15_Table)
CY15_Table[, which(colCY15_Table == "AGI")] <- CY_unionAGI #1
CY15_Table[match(TF_family$AGI, CY15_Table[, "AGI"]), which(colCY15_Table == "TF")] <- "Yes" #2

#######read_file#############################################################################################################################################################################
filename <- list.files("~/Nakano_RNAseq/network_analysis/base/change_MCLparameter/CY15", pattern=".csv", full.names = T)
temp2 <- rep(c("12h", "1h","24h", "3h"), each = 5)
temp3 <- rep(c("Botrytis_cinerea", "BTH", "MeJA", "PstDC3000", "sub1"), times = 4)
object_name <- c()
object_name_all <- c()
test <- c()
n <- 1
for (n in 1:length(filename)){
  object_name <- c(object_name, paste0("CY15", "_", temp2[n], "_", temp3[n]))
  assign(object_name[n], read.table(filename[n],  header=T, sep=",", stringsAsFactors = F))
  object_name_all <- c(object_name_all, object_name[n])
  n <- n+1
}

#####sub1_genelist##########################################
root <- list(CY15_12h_Botrytis_cinerea, CY15_12h_BTH, CY15_12h_MeJA, CY15_12h_PstDC3000,
             CY15_1h_Botrytis_cinerea, CY15_1h_BTH, CY15_1h_MeJA, CY15_1h_PstDC3000,
             CY15_24h_Botrytis_cinerea, CY15_24h_BTH, CY15_24h_MeJA, CY15_24h_PstDC3000,
             CY15_3h_Botrytis_cinerea, CY15_3h_BTH, CY15_3h_MeJA, CY15_3h_PstDC3000
)

root_name <- c("CY15_12h_Botrytis_cinerea", "CY15_12h_BTH", "CY15_12h_MeJA", "CY15_12h_PstDC3000",
               "CY15_1h_Botrytis_cinerea", "CY15_1h_BTH", "CY15_1h_MeJA", "CY15_1h_PstDC3000",
               "CY15_24h_Botrytis_cinerea", "CY15_24h_BTH", "CY15_24h_MeJA", "CY15_24h_PstDC3000",
               "CY15_3h_Botrytis_cinerea", "CY15_3h_BTH", "CY15_3h_MeJA", "CY15_3h_PstDC3000"
)
degree <- paste0(root_name, "_degree")
BC <- paste0(root_name, "_BetweennessCentrality")

root_name <- list(degree, BC)

rootob <- c("Degree", "BetweennessCentrality")
m <- 1
for(m in m:length(root_name)){
  n <- 1
  for(n in n:length(root_name[[m]])){
    temp <- match(root[[n]][, "name"], CY15_Table[, "AGI"])
    CY15_Table[temp, root_name[[m]][n]] <- root[[n]][, rootob[m]]
    n <- n+1
  }
  m <- m+1
}

root_name <- c("CY15_12h_Botrytis_cinerea_DEGs", "CY15_12h_BTH_DEGs", "CY15_12h_MeJA_DEGs", "CY15_12h_PstDC3000_DEGs",
               "CY15_1h_Botrytis_cinerea_DEGs", "CY15_1h_BTH_DEGs", "CY15_1h_MeJA_DEGs", "CY15_1h_PstDC3000_DEGs",
               "CY15_24h_Botrytis_cinerea_DEGs", "CY15_24h_BTH_DEGs", "CY15_24h_MeJA_DEGs", "CY15_24h_PstDC3000_DEGs",
               "CY15_3h_Botrytis_cinerea_DEGs", "CY15_3h_BTH_DEGs", "CY15_3h_MeJA_DEGs", "CY15_3h_PstDC3000_DEGs"
)

m <- 1
for(m in m:length(root_name)){
  temp <- match(root[[m]][, "name"], CY15_Table[, "AGI"])
  CY15_Table[temp, which(colCY15_Table == root_name[m])] <- "Yes"
  m <- m+1
}

#####################################################################################
root <- list(CY15_1h_sub1, CY15_3h_sub1, CY15_12h_sub1, CY15_24h_sub1)
root_name <- c("CY15_1h_sub1", "CY15_3h_sub1", "CY15_12h_sub1", "CY15_24h_sub1")
degree <- paste0(root_name, "_degree")
BC <- paste0(root_name, "_BetweennessCentrality")
root_name <- list(degree, BC)

rootob <- c("Degree", "BetweennessCentrality")
m <- 1
for(m in m:length(root_name)){
  n <- 1
  for(n in n:length(root_name[[m]])){
    temp <- match(root[[n]][, "name"], CY15_Table[, "AGI"])
    CY15_Table[temp, root_name[[m]][n]] <- root[[n]][, rootob[m]]
    n <- n+1
  }
  m <- m+1
}

root_name <- c("CY15_1h_sub1_DEGs", "CY15_3h_sub1_DEGs", "CY15_12h_sub1_DEGs", "CY15_24h_sub1_DEGs")
m <- 1
for(m in m:length(root_name)){
  temp <- match(root[[m]][, "name"], CY15_Table[, "AGI"])
  CY15_Table[temp, which(colCY15_Table == root_name[m])] <- "Yes"
  m <- m+1
}
#######check#############################################################################

library(dplyr)
cname <- c()
rank <- c()
rank_all <- list()
object <- c("degree", "BetweennessCentrality")
m <- 1
for(m in m:length(object)){
  temp <- CY15_Table[, grep(object[m], colnames(CY15_Table))]
  rownames(temp) <- CY15_Table[, "AGI"]
  cname <- c(cname, colnames(temp))
  n <- 1
  for(n in n:ncol(temp)){
    a <- as.vector(temp[, n])
    b <- rep(0, times = length(a))
    c <- which(a != "No")
    b[c] <- dense_rank(a[c])
    rank <- c(rank, b)
    n <- n+1
  }
  rank_all <- c(rank_all, list(rank))
  rank <- c()
  m <- m+1
}

ranking <- matrix(unlist(rank_all), nrow = length(b), ncol = length(object)*ncol(temp))
ranking[ranking == 0] <- ""
colnames(ranking) <- paste0(cname, "_ranking")
CY15_Table[CY15_Table == "No"] <- ""

CY15_Table <- data.frame(CY15_Table,
                         ranking)


write.table(CY15_Table, file = "~/Nakano_RNAseq/network_analysis/eachCY_Table/change_MCLparameter/CY15_Table_changeMCL.txt", append=F, quote = F, sep = "\t", row.names = F)