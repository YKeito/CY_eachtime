network_density <- read.table("~/Nakano_RNAseq/network_analysis/base/eachCY_131224h/network density/network density.txt", sep = "\t", header = T)
network_density <- data.frame(network_density, sort = rep(c(1:4), times =3))
colnames(network_density) <- c("CY", "value", "time", "sort")

g <- ggplot(network_density,
            aes(x = reorder(time, sort),
                y = value,
                group = CY
            )
)
g <- g + geom_line(aes(colour=CY))
g <- g + geom_point(aes(colour=CY))
g <- g + xlab("")
g <- g + ylab("network density")

plot(g)
ggsave(file = "Nakano_RNAseq/network_analysis/results_Fig/network_degree/network density.png", plot = g)